import argparse
import pathlib
import os
import time
import torch
import torchvision.utils as vutils
import numpy as np
import torchvision.models as models
from torchvision import datasets
from tensorboardX import SummaryWriter
from modules.python.Options import ImageSizeOptions
from modules.python.models.ModelHander import ModelHandler
from modules.python.models.simple_model import TransducerGRU
from modules.python.Options import ImageSizeOptions, TrainOptions
import torchvision


def visualize_model(model_path, output_path):
    writer = SummaryWriter(output_path)

    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    seq_len=ImageSizeOptions.SEQ_LENGTH,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS)
    dummy_input = torch.Tensor(1, TrainOptions.TRAIN_WINDOW, ImageSizeOptions.IMAGE_HEIGHT)
    hidden = torch.zeros(dummy_input.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
    # output = transducer_model(dummy_input, hidden)
    # exit()
    writer.add_graph(transducer_model, (dummy_input, hidden,))

    # dummy_input = torch.Tensor(1, 3, 224, 224)
    # model = torchvision.models.densenet121()
    # writer.add_graph(model, (dummy_input, ))
    writer.close()



def handle_output_directory(output_dir):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    timestr = time.strftime("%m%d%Y_%H%M%S")

    if output_dir[-1] != "/":
        output_dir += "/"

    output_dir += "tensorboard_output_"+str(timestr) + "/"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    return output_dir


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model_path",
        type=str,
        required=False,
        default='./model',
        help="Path of the model to load and retrain"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        default='./outputs/tensorboard/',
        help="Path and file_name to save model, default is ./model"
    )
    FLAGS, not_parsed = parser.parse_known_args()
    output_dir = handle_output_directory(FLAGS.output_dir)
    visualize_model(FLAGS.model_path, output_dir)
