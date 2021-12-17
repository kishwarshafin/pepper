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
from pepper_variant.modules.python.Options import ImageSizeOptions, ImageSizeOptionsHP, TrainOptions
from pepper_variant.modules.python.DataStorePredict import DataStore

os.environ['PYTHONWARNINGS'] = 'ignore:semaphore_tracker:UserWarning'


def predict(options, input_filepath, output_filepath, threads):
    if options.use_hp_info:
        image_features = ImageSizeOptionsHP.IMAGE_HEIGHT
    else:
        image_features = ImageSizeOptions.IMAGE_HEIGHT

    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(options.model_path,
                                                    image_features=image_features,
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
                             batch_size=options.batch_size,
                             shuffle=False,
                             num_workers=options.num_workers,
                             collate_fn=SequenceDataset.my_collate)
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
            outputs = transducer_model(images, hidden, cell_state, False)

            output_base, output_type = tuple(outputs)
            # output_base = output_base.detach().cpu().numpy()
            output_type = output_type.detach().cpu().numpy()

            prediction_data_file.write_prediction(batch_completed, contigs, positions, depths, candidates, candidate_frequencies, output_type)
            batch_completed += 1

            if batch_completed % 100 == 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " +
                                 "INFO: BATCHES PROCESSED " + str(batch_completed) + "/" + str(total_batches) + ".\n")
                sys.stderr.flush()


def predict_distributed_gpu(options, filepath, output_filepath, threads_per_caller):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param options: Options for prediction.
    :param filepath: Path to image files to predict on.
    :param output_filepath: Path to output directory.
    :param threads_per_caller: How many threads to set per caller.
    :return:
    """
    predict(options, filepath, output_filepath, threads_per_caller)
