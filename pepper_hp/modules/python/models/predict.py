import argparse
import sys
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from modules.python.models.dataloader_predict import SequenceDataset
from modules.python.TextColor import TextColor
from tqdm import tqdm
import numpy as np
from modules.python.models.ModelHander import ModelHandler
from modules.python.Options import ImageSizeOptions, TrainOptions
from modules.python.DataStorePredict import DataStore


def predict(test_file, output_filename, model_path, batch_size, threads, num_workers, gpu_mode):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param test_file: File to predict on
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param gpu_mode: If true, predictions will be done over GPU
    :param threads: Number of threads to set for pytorch
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
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
    sys.stderr.write(TextColor.CYAN + 'MODEL LOADED\n')

    with torch.no_grad():
        for contig, contig_start, contig_end, chunk_id, images, position, index, ref_seq, labels in tqdm(test_loader, ncols=50):
            sys.stderr.flush()
            images = images.type(torch.FloatTensor)

            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            prediction_base_counter = torch.zeros((images.size(0),
                                                   ImageSizeOptions.SEQ_LENGTH,
                                                   ImageSizeOptions.TOTAL_LABELS))

            if gpu_mode:
                images = images.cuda()
                hidden = hidden.cuda()
                prediction_base_counter = prediction_base_counter.cuda()

            for i in range(0, ImageSizeOptions.SEQ_LENGTH, TrainOptions.WINDOW_JUMP):
                if i + TrainOptions.TRAIN_WINDOW > ImageSizeOptions.SEQ_LENGTH:
                    break
                chunk_start = i
                chunk_end = i + TrainOptions.TRAIN_WINDOW

                # chunk all the data
                image_chunk = images[:, i:i+TrainOptions.TRAIN_WINDOW]

                # run inference
                out, hidden = transducer_model(image_chunk, hidden)

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
                    base_prediction = (inference_layers(out) * 10).type(torch.IntTensor).cuda()
                else:
                    base_prediction = (inference_layers(out) * 10).type(torch.IntTensor)

                # now simply add the tensor to the global counter
                prediction_base_counter = torch.add(prediction_base_counter, base_prediction)

            base_values, base_labels = torch.max(prediction_base_counter, 2)

            predicted_base_labels = base_labels.cpu().numpy()

            for i in range(images.size(0)):
                prediction_data_file.write_prediction(contig[i],
                                                      contig_start[i],
                                                      contig_end[i],
                                                      chunk_id[i],
                                                      position[i],
                                                      index[i],
                                                      ref_seq[i],
                                                      predicted_base_labels[i])
