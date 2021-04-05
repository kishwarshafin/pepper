import sys
import os
import torch
import time
import torch.onnx
import torch.nn as nn
import onnxruntime
from datetime import datetime
import concurrent.futures
from torch.utils.data import DataLoader
import numpy as np

from pepper_variant.modules.python.models.dataloader_predict import SequenceDatasetHP
from pepper_variant.modules.python.models.ModelHander import ModelHandler
from pepper_variant.modules.python.Options import ImageSizeOptionsHP, TrainOptions
from pepper_variant.modules.python.DataStorePredict import DataStore


def predict_hp(input_filepath, file_chunks, output_filepath, model_path, batch_size, num_workers, threads, thread_id):
    # session options
    sess_options = onnxruntime.SessionOptions()
    sess_options.intra_op_num_threads = threads
    sess_options.execution_mode = onnxruntime.ExecutionMode.ORT_SEQUENTIAL
    sess_options.graph_optimization_level = onnxruntime.GraphOptimizationLevel.ORT_ENABLE_ALL

    ort_session = onnxruntime.InferenceSession(model_path + ".onnx", sess_options=sess_options)
    torch.set_num_threads(threads)

    # create output file
    output_filename = output_filepath + "pepper_prediction_" + str(thread_id) + ".hdf"
    prediction_data_file = DataStore(output_filename, mode='w')

    # data loader
    input_data = SequenceDatasetHP(input_filepath, file_chunks)

    data_loader = DataLoader(input_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers)

    batch_completed = 0
    total_batches = len(data_loader)
    with torch.no_grad():
        for contig, contig_start, contig_end, chunk_id, images_hp1, images_hp2, position, index in data_loader:
            sys.stderr.flush()
            images_hp1 = images_hp1.type(torch.FloatTensor)
            images_hp2 = images_hp2.type(torch.FloatTensor)
            hidden_hp1 = torch.zeros(images_hp1.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
            hidden_hp2 = torch.zeros(images_hp2.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            prediction_base_tensor_hp1 = torch.zeros((images_hp1.size(0), images_hp1.size(1), ImageSizeOptionsHP.TOTAL_LABELS))
            prediction_base_tensor_hp2 = torch.zeros((images_hp2.size(0), images_hp2.size(1), ImageSizeOptionsHP.TOTAL_LABELS))

            for i in range(0, ImageSizeOptionsHP.SEQ_LENGTH, TrainOptions.WINDOW_JUMP):
                if i + TrainOptions.TRAIN_WINDOW > ImageSizeOptionsHP.SEQ_LENGTH:
                    break
                chunk_start = i
                chunk_end = i + TrainOptions.TRAIN_WINDOW
                # chunk all the data
                image_chunk_hp1 = images_hp1[:, chunk_start:chunk_end]
                image_chunk_hp2 = images_hp2[:, chunk_start:chunk_end]

                # run inference on onnx mode, which takes numpy inputs
                ort_inputs_hp1 = {ort_session.get_inputs()[0].name: image_chunk_hp1.cpu().numpy(),
                                  ort_session.get_inputs()[1].name: hidden_hp1.cpu().numpy()}
                ort_inputs_hp2 = {ort_session.get_inputs()[0].name: image_chunk_hp2.cpu().numpy(),
                                  ort_session.get_inputs()[1].name: hidden_hp2.cpu().numpy()}
                output_base_hp1, hidden_hp1 = ort_session.run(None, ort_inputs_hp1)
                output_base_hp2, hidden_hp2 = ort_session.run(None, ort_inputs_hp2)
                output_base_hp1 = torch.from_numpy(output_base_hp1)
                output_base_hp2 = torch.from_numpy(output_base_hp2)
                hidden_hp1 = torch.from_numpy(hidden_hp1)
                hidden_hp2 = torch.from_numpy(hidden_hp2)

                # now calculate how much padding is on the top and bottom of this chunk so we can do a simple
                # add operation
                top_zeros = chunk_start
                bottom_zeros = ImageSizeOptionsHP.SEQ_LENGTH - chunk_end

                # do softmax and get prediction
                # we run a softmax a padding to make the output tensor compatible for adding
                inference_layers = nn.Sequential(
                    nn.Softmax(dim=2),
                    nn.ZeroPad2d((0, 0, top_zeros, bottom_zeros))
                )

                # run the softmax and padding layers
                base_prediction_hp1 = (inference_layers(output_base_hp1) * 10000).type(torch.IntTensor)
                base_prediction_hp2 = (inference_layers(output_base_hp2) * 10000).type(torch.IntTensor)

                # now simply add the tensor to the global counter
                prediction_base_tensor_hp1 = torch.add(prediction_base_tensor_hp1, base_prediction_hp1)
                prediction_base_tensor_hp2 = torch.add(prediction_base_tensor_hp2, base_prediction_hp2)

            # base_values, base_labels = torch.max(prediction_base_tensor, 2)
            #
            # predicted_base_labels = base_labels.cpu().numpy()
            prediction_base_tensor_hp1 = prediction_base_tensor_hp1.cpu().numpy().astype(int)
            prediction_base_tensor_hp2 = prediction_base_tensor_hp2.cpu().numpy().astype(int)

            for i in range(images_hp1.size(0)):
                prediction_data_file.write_prediction_hp(contig[i],
                                                         contig_start[i],
                                                         contig_end[i],
                                                         chunk_id[i],
                                                         position[i],
                                                         index[i],
                                                         prediction_base_tensor_hp1[i],
                                                         prediction_base_tensor_hp2[i])
            batch_completed += 1

            if thread_id == 0 and batch_completed % 5 == 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " +
                                 "INFO: BATCHES PROCESSED " + str(batch_completed) + "/" + str(total_batches) + ".\n")
                sys.stderr.flush()

    return thread_id


def predict_hp_distributed_cpu(filepath, file_chunks, output_filepath, model_path, batch_size, total_callers, threads, num_workers):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param filepath: Path to image files to predict on
    :param file_chunks: Path to chunked files
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param output_filepath: Path to output directory
    :param total_callers: Number of callers to spawn
    :param threads: Number of threads to use per caller
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    # load the model and create an ONNX session
    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    input_channels=ImageSizeOptionsHP.IMAGE_CHANNELS,
                                                    image_features=ImageSizeOptionsHP.IMAGE_HEIGHT,
                                                    seq_len=ImageSizeOptionsHP.SEQ_LENGTH,
                                                    num_classes=ImageSizeOptionsHP.TOTAL_LABELS)
    transducer_model.eval()

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL LOADING TO ONNX\n")
    x = torch.zeros(1, TrainOptions.TRAIN_WINDOW, ImageSizeOptionsHP.IMAGE_HEIGHT)
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
        futures = [executor.submit(predict_hp, filepath, file_chunks[thread_id], output_filepath, model_path, batch_size, num_workers, threads, thread_id)
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


