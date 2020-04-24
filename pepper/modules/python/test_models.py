import os
import sys

import torch
import torch.nn.parallel
from datetime import datetime
from pepper.modules.python.models.test import test
from pepper.modules.python.models.ModelHander import ModelHandler
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
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]  Loading data\n")

    if os.path.isfile(model_path) is False:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: INVALID PATH TO MODEL\n")
        exit(1)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]  INFO: MODEL LOADING\n")

    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    seq_len=ImageSizeOptions.SEQ_LENGTH,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]  INFO: MODEL LOADED\n")
    sys.stderr.flush()

    if gpu_mode:
        transducer_model = torch.nn.DataParallel(transducer_model).cuda()

    stats_dictioanry = test(test_file, batch_size, gpu_mode, transducer_model, num_workers,
                            gru_layers, hidden_size, num_classes=ImageSizeOptions.TOTAL_LABELS,
                            print_details=print_details)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST COMPLETE")


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
