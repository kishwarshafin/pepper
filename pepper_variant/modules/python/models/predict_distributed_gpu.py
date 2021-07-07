import sys
import os
import torch
import torch.distributed as dist
import torch.nn as nn
from torch.utils.data import DataLoader
import torch.multiprocessing as mp
from torch.nn.parallel import DataParallel
from datetime import datetime

from pepper_variant.modules.python.models.dataloader_predict import SequenceDataset
from pepper_variant.modules.python.models.ModelHander import ModelHandler
from pepper_variant.modules.python.Options import ImageSizeOptions, TrainOptions
from pepper_variant.modules.python.DataStorePredict import DataStore

os.environ['PYTHONWARNINGS'] = 'ignore:semaphore_tracker:UserWarning'


def predict(input_filepath, output_filepath, model_path, batch_size, num_workers, threads):
    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS,
                                                    num_type_classes=ImageSizeOptions.TOTAL_TYPE_LABELS)

    transducer_model.eval()
    transducer_model = transducer_model.eval()
    # create output file
    output_filename = output_filepath + "pepper_prediction" + ".hdf"
    prediction_data_file = DataStore(output_filename, mode='w')

    # data loader
    input_data = SequenceDataset(input_filepath)
    data_loader = DataLoader(input_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers)
    torch.set_num_threads(threads)

    transducer_model = transducer_model.cuda()
    transducer_model.eval()
    transducer_model = DataParallel(transducer_model)

    batch_completed = 0
    total_batches = len(data_loader)

    with torch.no_grad():
        for contigs, positions, depths, candidates, candidate_frequencies, images in data_loader:
            sys.stderr.flush()
            images = images.type(torch.FloatTensor)
            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
            cell_state = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            images = images.cuda()
            hidden = hidden.cuda()
            cell_state = cell_state.cuda()

            # run inference
            # output_base, output_type = transducer_model(images, hidden, cell_state, False)
            output_base = transducer_model(images, hidden, cell_state, False)

            output_base = output_base.detach().cpu().numpy()
            # output_type = output_type.detach().cpu().numpy()

            prediction_data_file.write_prediction(batch_completed, contigs, positions, depths, candidates, candidate_frequencies, output_base)
            batch_completed += 1

            if batch_completed % 5 == 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " +
                                 "INFO: BATCHES PROCESSED " + str(batch_completed) + "/" + str(total_batches) + ".\n")
                sys.stderr.flush()


# def cleanup():
#     dist.destroy_process_group()
#
#
# def setup(rank, total_callers, args, all_input_files):
#     os.environ['MASTER_ADDR'] = 'localhost'
#     os.environ['MASTER_PORT'] = '12355'
#
#     # initialize the process group
#     dist.init_process_group("gloo", rank=rank, world_size=total_callers)
#
#     filepath, output_filepath, model_path, batch_size, threads_per_caller, device_ids, num_workers = args
#
#     # issue with semaphore lock: https://github.com/pytorch/pytorch/issues/2517
#     # mp.set_start_method('spawn')
#
#     # Explicitly setting seed to make sure that models created in two processes
#     # start from same random weights and biases. https://github.com/pytorch/pytorch/issues/2517
#     # torch.manual_seed(42)
#     predict(filepath, all_input_files[rank],  output_filepath, model_path, batch_size, num_workers, threads_per_caller, device_ids[rank], rank)
#     cleanup()


def predict_distributed_gpu(filepath, file_chunks, output_filepath, model_path, batch_size, total_callers, threads_per_caller, device_ids, num_workers):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param filepath: Path to image files to predict on
    :param file_chunks: Path to chunked files
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param output_filepath: Path to output directory
    :param total_callers: Number of callers
    :param threads_per_caller: How many threads to set per caller
    :param device_ids: Device ID of GPU to be used
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    # args = (filepath, output_filepath, model_path, batch_size, threads_per_caller, device_ids, num_workers)

    predict(filepath, output_filepath, model_path, batch_size, num_workers, threads_per_caller)
    # mp.spawn(setup,
    #          args=(total_callers, args, file_chunks),
    #          nprocs=total_callers,
    #          join=True)
