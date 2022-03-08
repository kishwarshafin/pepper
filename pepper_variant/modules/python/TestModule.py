import os
import sys
import torch
import torch.nn.parallel
from datetime import datetime
from pepper_variant.modules.python.models.test import test
from pepper_variant.modules.python.models.test_hp import test_hp
from pepper_variant.modules.python.models.ModelHander import ModelHandler
from pepper_variant.modules.python.Options import ImageSizeOptions, ImageSizeOptionsHP


def do_test(test_file, use_hp_info, batch_size, gpu_mode, num_workers, model_path):
    """
    Train a model and save
    :param test_file: A CSV file containing test image information
    :param batch_size: Batch size for training
    :param gpu_mode: If true the model will be trained on GPU
    :param num_workers: Number of workers for data loading
    :param model_path: Path to a saved model
    :param num_classes: Number of output classes
    :return:
    """
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: LOADING DATA\n")

    if os.path.isfile(model_path) is False:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: INVALID PATH TO MODEL\n")
        exit(1)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL LOADING\n")

    if not use_hp_info:
        image_height = ImageSizeOptions.IMAGE_HEIGHT
        num_classes = ImageSizeOptions.TOTAL_LABELS
        num_type_classes = ImageSizeOptions.TOTAL_TYPE_LABELS
    else:
        image_height = ImageSizeOptionsHP.IMAGE_HEIGHT
        num_classes = ImageSizeOptionsHP.TOTAL_LABELS
        num_type_classes = ImageSizeOptionsHP.TOTAL_TYPE_LABELS

    if use_hp_info:
        transducer_model, hidden_size, gru_layers, prev_ite = \
            ModelHandler.load_simple_model_for_training(model_path,
                                                        image_features=image_height,
                                                        num_classes=num_classes,
                                                        num_type_classes=num_type_classes)
    else:
        transducer_model, hidden_size, gru_layers, prev_ite = \
            ModelHandler.load_simple_model_for_training(model_path,
                                                        image_features=image_height,
                                                        num_classes=num_classes,
                                                        num_type_classes=num_type_classes)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL LOADED\n")

    if gpu_mode:
        transducer_model = torch.nn.DataParallel(transducer_model).cuda()

    if use_hp_info:
        stats_dictionary = test(test_file, batch_size, gpu_mode, transducer_model, num_workers)
    else:
        stats_dictionary = test(test_file, batch_size, gpu_mode, transducer_model, num_workers)