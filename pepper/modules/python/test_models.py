import os
import sys

import torch
import torch.nn.parallel
from pepper.modules.python.models.test import test
from pepper.modules.python.models.ModelHander import ModelHandler
from pepper.modules.python.TextColor import TextColor
from pepper.modules.python.Options import ImageSizeOptions
"""
TEST A TRAINED MODEL
"""


def do_test(test_file, batch_size, gpu_mode, num_workers, model_path, print_details):
    """
    Train a model and save
    :param test_file: A CSV file containing test image information
    :param batch_size: Batch size for training
    :param gpu_mode: If true the model will be trained on GPU
    :param num_workers: Number of workers for data loading
    :param model_path: Path to a saved model
    :param print_details: Print debug stuff
    :return:
    """
    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)

    if os.path.isfile(model_path) is False:
        sys.stderr.write(TextColor.RED + "ERROR: INVALID PATH TO MODEL\n")
        exit(1)

    sys.stderr.write(TextColor.GREEN + "INFO: MODEL LOADING\n" + TextColor.END)

    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    seq_len=ImageSizeOptions.SEQ_LENGTH,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS)

    sys.stderr.write(TextColor.GREEN + "INFO: MODEL LOADED\n" + TextColor.END)

    if gpu_mode:
        transducer_model = torch.nn.DataParallel(transducer_model).cuda()

    stats_dictioanry = test(test_file, batch_size, gpu_mode, transducer_model, num_workers,
                            gru_layers, hidden_size, num_classes=ImageSizeOptions.TOTAL_LABELS,
                            print_details=print_details)

    sys.stderr.write(TextColor.PURPLE + 'DONE\n' + TextColor.END)


def test_models(image_dir, model_path, num_workers, batch_size, gpu_mode, print_details):
    """
    Processes arguments and performs tasks.
    """
    do_test(image_dir,
            batch_size,
            gpu_mode,
            num_workers,
            model_path,
            print_details)
