import os
import time
from datetime import datetime
import torch
import sys

from pepper_snp.modules.python.models.train import train
from pepper_snp.modules.python.models.train_distributed import train_distributed
from pepper_snp.modules.python.Options import TrainOptions


class TrainModule:
    """
    Train module
    """
    def __init__(self, train_file, test_file, gpu_mode, max_epochs, batch_size, num_workers,
                 retrain_model, retrain_model_path, model_dir, stats_dir):
        self.train_file = train_file
        self.test_file = test_file
        self.gpu_mode = gpu_mode
        self.log_directory = stats_dir
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

    def train_model_distributed(self, device_ids, callers_per_gpu):
        """
        DO DISTRIBUTED GPU INFERENCE. THIS MODE WILL ENABLE ONE MODEL PER GPU
        """
        start_time = time.time()
        if not torch.cuda.is_available():
            sys.stderr.write("ERROR: TORCH IS NOT BUILT WITH CUDA.\n")
            sys.stderr.write("SEE TORCH CAPABILITY:\n$ python3\n  "
                             ">>> import torch \n "
                             ">>> torch.cuda.is_available()\n "
                             "If true then cuda is avilable")
            exit(1)

        if device_ids is None:
            total_gpu_devices = torch.cuda.device_count()
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL GPU AVAILABLE: " + str(total_gpu_devices) + "\n")
            device_ids = [i for i in range(0, total_gpu_devices)]
        else:
            device_ids = [int(i) for i in device_ids.split(',')]
            device_ids = list(set(device_ids))
            for device_id in device_ids:
                major_capable, minor_capable = torch.cuda.get_device_capability(device=device_id)
                if major_capable < 0:
                    sys.stderr.write("ERROR: GPU DEVICE: " + str(device_id) + " IS NOT CUDA CAPABLE.\n")
                    sys.stderr.write("Try running: $ python3\n"
                                     ">>> import torch \n"
                                     ">>> torch.cuda.get_device_capability(device=" + str(device_id) + ")\n")
                    exit(1)
                else:
                    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CAPABILITY OF GPU#"
                                     + str(device_id) + ":\t" + str(major_capable)
                                     + "-" + str(minor_capable) + "\n")

        multiplied_device_ids = []
        for device_id in device_ids:
            for i in range(0, callers_per_gpu):
                multiplied_device_ids.append(device_id)
        device_ids = multiplied_device_ids

        total_callers = len(device_ids)

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: AVAILABLE GPU DEVICES: " + str(list(set(device_ids))) + "\n")
        if total_callers == 0:
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: NO GPU AVAILABLE BUT GPU MODE IS SET\n")
            exit(1)
            total_gpu_devices = torch.cuda.device_count()

            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL GPU AVAILABLE: " + str(total_gpu_devices) + "\n")
            device_ids = [i for i in range(0, total_gpu_devices)]
            total_callers = total_gpu_devices

            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: AVAILABLE GPU DEVICES: " + str(device_ids) + "\n")

            if total_callers == 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: NO GPU AVAILABLE BUT GPU MODE IS SET\n")
                exit()

        # train a model
        train_distributed(self.train_file,
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
                          device_ids,
                          total_callers,
                          train_mode=True)


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


def train_pepper_snp_model(train_file, test_file, output_dir, gpu_mode, epoch_size, batch_size, num_workers,
                           retrain_model, retrain_model_path, distributed, device_ids, per_gpu_callers):
    model_out_dir, log_dir = handle_output_directory(output_dir)
    tm = TrainModule(train_file,
                     test_file,
                     gpu_mode,
                     epoch_size,
                     batch_size,
                     num_workers,
                     retrain_model,
                     retrain_model_path,
                     model_out_dir,
                     log_dir)

    if distributed:
        tm.train_model_distributed(device_ids, per_gpu_callers)
    else:
        tm.train_model()