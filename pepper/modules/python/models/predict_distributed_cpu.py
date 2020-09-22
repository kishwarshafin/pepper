import onnxruntime
import sys
import os
import time
import torch
import torch.nn as nn
import concurrent.futures
from datetime import datetime
from torch.utils.data import DataLoader
from pepper.modules.python.models.dataloader_predict import SequenceDataset
from pepper.modules.python.models.ModelHander import ModelHandler
from pepper.modules.python.Options import ImageSizeOptions, TrainOptions
from pepper.modules.python.DataStorePredict import DataStore
import torch.onnx
import warnings
warnings.filterwarnings("ignore", message="Exporting a model to ONNX with a batch_size other than 1, with a variable lenght with GRU can cause an error when running the ONNX model with a different batch size")


def predict(input_filepath, file_chunks, output_filepath, batch_size, num_workers, rank, threads_per_caller, model_path):
    # session options
    sess_options = onnxruntime.SessionOptions()
    sess_options.intra_op_num_threads = threads_per_caller
    sess_options.execution_mode = onnxruntime.ExecutionMode.ORT_SEQUENTIAL
    sess_options.graph_optimization_level = onnxruntime.GraphOptimizationLevel.ORT_ENABLE_ALL

    ort_session = onnxruntime.InferenceSession(model_path + ".onnx", sess_options=sess_options)
    torch.set_num_threads(threads_per_caller)

    # create output file
    output_filename = output_filepath + "pepper_prediction_" + str(rank) + ".hdf"
    prediction_data_file = DataStore(output_filename, mode='w')

    # data loader
    input_data = SequenceDataset(input_filepath, file_chunks)
    data_loader = DataLoader(input_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers)

    batch_completed = 0

    with torch.no_grad():
        for contig, contig_start, contig_end, chunk_id, images, position, index in data_loader:
            images = images.type(torch.FloatTensor)

            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            prediction_base_tensor = torch.zeros((images.size(0), images.size(1), ImageSizeOptions.TOTAL_LABELS))

            for i in range(0, ImageSizeOptions.SEQ_LENGTH, TrainOptions.WINDOW_JUMP):
                if i + TrainOptions.TRAIN_WINDOW > ImageSizeOptions.SEQ_LENGTH:
                    break
                chunk_start = i
                chunk_end = i + TrainOptions.TRAIN_WINDOW
                # chunk all the data
                image_chunk = images[:, chunk_start:chunk_end]

                # run inference on onnx mode, which takes numpy inputs
                ort_inputs = {ort_session.get_inputs()[0].name: image_chunk.cpu().numpy(),
                              ort_session.get_inputs()[1].name: hidden.cpu().numpy()}
                output_base, hidden = ort_session.run(None, ort_inputs)
                output_base = torch.from_numpy(output_base)
                hidden = torch.from_numpy(hidden)

                # now calculate how much padding is on the top and bottom of this chunk so we can do a simple
                # add operation
                top_zeros = chunk_start
                bottom_zeros = ImageSizeOptions.SEQ_LENGTH - chunk_end

                counts = torch.ones((output_base.size(0), output_base.size(1), 1))

                # do softmax and get prediction
                # we run a softmax a padding to make the output tensor compatible for adding
                inference_layers = nn.Sequential(
                    nn.Softmax(dim=2),
                    nn.ZeroPad2d((0, 0, top_zeros, bottom_zeros))
                )
                base_prediction = inference_layers(output_base)

                # now simply add the tensor to the global counter
                prediction_base_tensor = torch.add(prediction_base_tensor, base_prediction)

            base_values, base_labels = torch.max(prediction_base_tensor, 2)

            # this part is for the phred score calculation
            counts = torch.ones((base_values.size(0), base_values.size(1) - 2 * ImageSizeOptions.SEQ_OVERLAP))
            top_ones = nn.ZeroPad2d((ImageSizeOptions.SEQ_OVERLAP, ImageSizeOptions.SEQ_OVERLAP))
            counts = top_ones(counts) + 1
            phred_score = -10 * torch.log10(1.0 - (base_values / counts))
            phred_score[phred_score == float('inf')] = 100

            predicted_base_labels = base_labels.cpu().numpy()
            phred_score = phred_score.cpu().numpy()

            for i in range(images.size(0)):
                prediction_data_file.write_prediction(contig[i], contig_start[i], contig_end[i], chunk_id[i],
                                                      position[i], index[i], predicted_base_labels[i], phred_score[i])

            if rank == 0:
                batch_completed += 1
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " +
                                 "INFO: BATCHES PROCESSED " + str(batch_completed) + "/" + str(len(data_loader)) + ".\n")
                sys.stderr.flush()

    return rank


def predict_cpu(filepath, file_chunks, output_filepath, model_path, batch_size, total_callers, threads_per_caller, num_workers):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param filepath: Path to image files to predict on
    :param file_chunks: Path to chunked files
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param output_filepath: Path to output directory
    :param total_callers: Number of callers to spawn
    :param threads_per_caller: Number of threads to use per caller
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    # load the model and create an ONNX session
    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    seq_len=ImageSizeOptions.SEQ_LENGTH,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS)
    transducer_model.eval()

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL LOADING TO ONNX\n")
    x = torch.zeros(1, TrainOptions.TRAIN_WINDOW, ImageSizeOptions.IMAGE_HEIGHT)
    h = torch.zeros(1, 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

    if not os.path.isfile(model_path + ".onnx"):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: SAVING MODEL TO ONNX\n")
        torch.onnx.export(transducer_model, (x, h),
                          model_path + ".onnx",
                          training=False,
                          opset_version=10,
                          do_constant_folding=True,
                          input_names=['input_image', 'input_hidden'],
                          output_names=['output_pred', 'output_hidden'],
                          dynamic_axes={'input_image': {0: 'batch_size'},
                                        'input_hidden': {0: 'batch_size'},
                                        'output_pred': {0: 'batch_size'},
                                        'output_hidden': {0: 'batch_size'}})

    start_time = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=total_callers) as executor:
        futures = [executor.submit(predict, filepath, file_chunks[thread_id], output_filepath, batch_size, num_workers, thread_id, threads_per_caller, model_path)
                   for thread_id in range(0, total_callers)]

        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                # get the results
                thread_id = fut.result()
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: THREAD "
                                 + str(thread_id) + " FINISHED SUCCESSFULLY.\n")
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    end_time = time.time()
    mins = int((end_time - start_time) / 60)
    secs = int((end_time - start_time)) % 60
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PREDICTION\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec\n")
