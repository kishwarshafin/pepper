import argparse
import os
import time

# Custom generator for our dataset
from modules.python.models.train import train
from modules.python.Options import TrainOptions
"""
Input:
- A directory containing labeled HDF5 files
- A directory containing labeled HDF5 files

Output:
- A set of trained models, one per epoch
"""


class TrainModule:
    """
    Train module
    """
    def __init__(self, train_file, test_file, gpu_mode, max_epochs, batch_size, num_workers,
                 retrain_model, retrain_model_path, model_dir, log_dir, stats_dir):
        self.train_file = train_file
        self.test_file = test_file
        self.gpu_mode = gpu_mode
        self.log_directory = log_dir
        self.model_dir = model_dir
        self.epochs = max_epochs
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.retrain_model = retrain_model
        self.retrain_model_path = retrain_model_path
        self.stats_dir = stats_dir
        self.hidden_size = TrainOptions.HIDDEN_SIZE
        self.gru_layers = TrainOptions.GRU_LAYERS
        # {'l2': 1.4946789574136535e-05, 'lr': 0.000541365592065579}
        self.learning_rate = 0.0001
        self.weight_decay = 0

    def train_model(self):
        # train a model
        model, optimizer, stats_dictionary = train(self.train_file,
                                                   self.test_file,
                                                   self.batch_size,
                                                   self.epochs,
                                                   self.gpu_mode,
                                                   self.num_workers,
                                                   self.retrain_model,
                                                   self.retrain_model_path,
                                                   self.gru_layers,
                                                   self.hidden_size,
                                                   self.learning_rate,
                                                   self.weight_decay,
                                                   self.model_dir,
                                                   self.stats_dir,
                                                   train_mode=True)

        return model, optimizer, stats_dictionary


def handle_output_directory(output_dir):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    timestr = time.strftime("%m%d%Y_%H%M%S")
    # process the output directory
    if output_dir[-1] != "/":
        output_dir += "/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # create an internal directory so we don't overwrite previous runs
    model_save_dir = output_dir + "trained_models_" + timestr + "/"
    if not os.path.exists(model_save_dir):
        os.mkdir(model_save_dir)

    stats_directory = model_save_dir + "stats_" + timestr + "/"

    if not os.path.exists(stats_directory):
        os.mkdir(stats_directory)

    return model_save_dir, stats_directory


def train_models(train_image_dir,
                 test_image_dir,
                 output_directory,
                 epoch_size,
                 batch_size,
                 num_workers,
                 retrain_model,
                 retrain_model_path,
                 gpu_mode):

    model_out_dir, log_dir = handle_output_directory(output_directory)
    tm = TrainModule(train_image_dir,
                     test_image_dir,
                     gpu_mode,
                     epoch_size,
                     batch_size,
                     num_workers,
                     retrain_model,
                     retrain_model_path,
                     model_out_dir,
                     log_dir,
                     log_dir)
    tm.train_model()
