import sys
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from modules.python.models.dataloader_predict import SequenceDataset
from modules.python.TextColor import TextColor
from tqdm import tqdm
from modules.python.models.ModelHander import ModelHandler
from modules.python.Options import ImageSizeOptions, TrainOptions
from modules.python.DataStorePredict import DataStore
from torch.utils import mkldnn
import torch.onnx
import os
import onnx
import onnxruntime


def to_numpy(tensor):
    return tensor.detach().cpu().numpy() if tensor.requires_grad else tensor.cpu().numpy()


def predict(test_file, output_filename, model_path, batch_size, threads, num_workers, gpu_mode, onnx_mode):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param test_file: File to predict on
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param gpu_mode: If true, predictions will be done over GPU
    :param onnx_mode: If true, onnx mode is on for cpu
    :param threads: Number of threads to set for pytorch
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    if gpu_mode:
        onnx_mode = False

    prediction_data_file = DataStore(output_filename, mode='w')
    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)

    torch.set_num_threads(threads)
    sys.stderr.write(TextColor.GREEN + 'INFO: TORCH THREADS SET TO: ' + str(torch.get_num_threads()) + ".\n"
                     + TextColor.END)

    sys.stderr.flush()

    # data loader
    test_data = SequenceDataset(test_file)
    test_loader = DataLoader(test_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers)

    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    seq_len=ImageSizeOptions.SEQ_LENGTH,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS)
    transducer_model.eval()

    if gpu_mode:
        transducer_model = torch.nn.DataParallel(transducer_model).cuda()
    elif onnx_mode:
        sys.stderr.write("INFO: MODEL LOADING TO ONNX\n")
        x = torch.zeros(1, TrainOptions.TRAIN_WINDOW, ImageSizeOptions.IMAGE_HEIGHT)
        h = torch.zeros(1, 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

        if not os.path.isfile(model_path + ".onnx"):
            sys.stderr.write("INFO: SAVING MODEL TO ONNX\n")
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

        sys.stderr.write("INFO: LOADING ONNX MODEL\n")
        onnx_model = onnx.load(model_path + ".onnx")
        # sys.stderr.write("INFO: CHECKING ONNX MODEL\n")
        # onnx.checker.check_model(onnx_model)
        sys.stderr.write("INFO: ONNX SESSION INITIALIZING\n")
        ort_session = onnxruntime.InferenceSession(model_path + ".onnx")

    sys.stderr.write(TextColor.CYAN + 'STARTING INFERENCE\n' + TextColor.END)

    with torch.no_grad():
        for contig, contig_start, contig_end, chunk_id, images, position, index in tqdm(test_loader, ncols=50):
            images = images.type(torch.FloatTensor)

            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            prediction_base_tensor = torch.zeros((images.size(0), images.size(1), ImageSizeOptions.TOTAL_LABELS))

            if gpu_mode:
                images = images.cuda()
                hidden = hidden.cuda()
                prediction_base_tensor = prediction_base_tensor.cuda()

            for i in range(0, ImageSizeOptions.SEQ_LENGTH, TrainOptions.WINDOW_JUMP):
                if i + TrainOptions.TRAIN_WINDOW > ImageSizeOptions.SEQ_LENGTH:
                    break
                chunk_start = i
                chunk_end = i + TrainOptions.TRAIN_WINDOW
                # chunk all the data
                image_chunk = images[:, chunk_start:chunk_end]

                if onnx_mode:
                    # run inference on onnx mode, which takes numpy inputs
                    ort_inputs = {ort_session.get_inputs()[0].name: image_chunk.cpu().numpy(),
                                  ort_session.get_inputs()[1].name: hidden.cpu().numpy()}
                    output_base, hidden = ort_session.run(None, ort_inputs)
                    output_base = torch.from_numpy(output_base)
                    hidden = torch.from_numpy(hidden)
                else:
                    # run inference
                    output_base, hidden = transducer_model(image_chunk, hidden)

                # now calculate how much padding is on the top and bottom of this chunk so we can do a simple
                # add operation
                top_zeros = chunk_start
                bottom_zeros = ImageSizeOptions.SEQ_LENGTH - chunk_end

                # do softmax and get prediction
                # we run a softmax a padding to make the output tensor compatible for adding
                inference_layers = nn.Sequential(
                    nn.Softmax(dim=2),
                    nn.ZeroPad2d((0, 0, top_zeros, bottom_zeros))
                )
                if gpu_mode:
                    inference_layers = inference_layers.cuda()

                # run the softmax and padding layers
                if gpu_mode:
                    base_prediction = inference_layers(output_base).cuda()
                else:
                    base_prediction = inference_layers(output_base)

                # now simply add the tensor to the global counter
                prediction_base_tensor = torch.add(prediction_base_tensor, base_prediction)

            base_values, base_labels = torch.max(prediction_base_tensor, 2)

            predicted_base_labels = base_labels.cpu().numpy()

            for i in range(images.size(0)):
                prediction_data_file.write_prediction(contig[i], contig_start[i], contig_end[i], chunk_id[i],
                                                      position[i], index[i], predicted_base_labels[i])
