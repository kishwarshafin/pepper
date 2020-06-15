import sys
import torch
import os
from tqdm import tqdm
import torch.distributed as dist
import torch.nn as nn
import torch.multiprocessing as mp
# Custom generator for our dataset
from datetime import datetime
from torch.utils.data import DataLoader
from pepper.modules.python.models.dataloader import SequenceDataset
from pepper.modules.python.models.ModelHander import ModelHandler
from pepper.modules.python.models.test import test
from pepper.modules.python.Options import ImageSizeOptions, TrainOptions


os.environ['PYTHONWARNINGS'] = 'ignore:semaphore_tracker:UserWarning'

"""
Train a model and return the model and optimizer trained.

Input:
- A train CSV containing training image set information (usually chr1-18)

Return:
- A trained model
"""
CLASS_WEIGHTS = [1.0, 1.0, 1.0, 1.0, 1.0]


def save_best_model(transducer_model, model_optimizer, hidden_size, layers, epoch,
                    file_name):
    """
    Save the best model
    :param transducer_model: A trained model
    :param model_optimizer: Model optimizer
    :param hidden_size: Number of hidden layers
    :param layers: Number of GRU layers to use
    :param epoch: Epoch/iteration number
    :param file_name: Output file name
    :return:
    """
    if os.path.isfile(file_name):
        os.remove(file_name)
    ModelHandler.save_checkpoint({
        'model_state_dict': transducer_model.state_dict(),
        'model_optimizer': model_optimizer.state_dict(),
        'hidden_size': hidden_size,
        'gru_layers': layers,
        'epochs': epoch,
    }, file_name)
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL SAVED SUCCESSFULLY.\n")


def train(train_file, test_file, batch_size, epoch_limit, gpu_mode, num_workers, retrain_model,
          retrain_model_path, gru_layers, hidden_size, lr, decay, model_dir, stats_dir, train_mode,
          world_size, rank, device_id):

    if train_mode is True and rank == 0:
        train_loss_logger = open(stats_dir + "train_loss.csv", 'w')
        test_loss_logger = open(stats_dir + "test_loss.csv", 'w')
        confusion_matrix_logger = open(stats_dir + "confusion_matrix.txt", 'w')
    else:
        train_loss_logger = None
        test_loss_logger = None
        confusion_matrix_logger = None

    torch.cuda.set_device(device_id)

    if rank == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: LOADING DATA\n")

    train_data_set = SequenceDataset(train_file)

    train_sampler = torch.utils.data.distributed.DistributedSampler(
        train_data_set,
        num_replicas=world_size,
        rank=rank
    )

    train_loader = torch.utils.data.DataLoader(
        dataset=train_data_set,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0,
        pin_memory=True,
        sampler=train_sampler)

    num_classes = ImageSizeOptions.TOTAL_LABELS

    if retrain_model is True:
        if os.path.isfile(retrain_model_path) is False:
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: INVALID PATH TO RETRAIN PATH MODEL --retrain_model_path\n")
            exit(1)
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RETRAIN MODEL LOADING\n")
        transducer_model, hidden_size, gru_layers, prev_ite = \
            ModelHandler.load_simple_model_for_training(retrain_model_path,
                                                        input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                        image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                        seq_len=ImageSizeOptions.SEQ_LENGTH,
                                                        num_classes=num_classes)

        if train_mode is True:
            epoch_limit = prev_ite + epoch_limit

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RETRAIN MODEL LOADED\n")
    else:
        transducer_model = ModelHandler.get_new_gru_model(input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                          image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                          gru_layers=gru_layers,
                                                          hidden_size=hidden_size,
                                                          num_classes=num_classes)
        prev_ite = 0

    param_count = sum(p.numel() for p in transducer_model.parameters() if p.requires_grad)
    if rank == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL TRAINABLE PARAMETERS:\t" + str(param_count) + "\n")

    model_optimizer = torch.optim.Adam(transducer_model.parameters(), lr=lr, weight_decay=decay)

    if retrain_model is True:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: OPTIMIZER LOADING\n")
        model_optimizer = ModelHandler.load_simple_optimizer(model_optimizer, retrain_model_path, gpu_mode)
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: OPTIMIZER LOADED\n")

    if gpu_mode:
        transducer_model = transducer_model.to(device_id)
        transducer_model = nn.parallel.DistributedDataParallel(transducer_model, device_ids=[device_id])

    class_weights = torch.Tensor(CLASS_WEIGHTS)
    # Loss
    criterion = nn.CrossEntropyLoss(class_weights)

    if gpu_mode is True:
        criterion = criterion.to(device_id)

    start_epoch = prev_ite

    # Train the Model
    if rank == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: Training starting\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] Start: " + str(start_epoch + 1) + " End: " + str(epoch_limit) + "\n")

    stats = dict()
    stats['loss_epoch'] = []
    stats['accuracy_epoch'] = []

    for epoch in range(start_epoch, epoch_limit, 1):
        total_loss = 0
        total_images = 0
        if rank == 0:
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] Train epoch: " + str(epoch + 1) + "\n")
        # make sure the model is in train mode. BN is different in train and eval.

        batch_no = 1
        if rank == 0:
            progress_bar = tqdm(
                total=len(train_loader),
                ncols=100,
                leave=True,
                position=rank,
                desc="Loss: ",
            )
        else:
            progress_bar = None

        transducer_model.train()
        for images, labels in train_loader:
            labels = labels.type(torch.LongTensor)
            images = images.type(torch.FloatTensor)
            if gpu_mode:
                images = images.to(device_id)
                labels = labels.to(device_id)

            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            if gpu_mode:
                hidden = hidden.to(device_id)

            for i in range(0, ImageSizeOptions.SEQ_LENGTH, TrainOptions.WINDOW_JUMP):
                model_optimizer.zero_grad()

                if i + TrainOptions.TRAIN_WINDOW > ImageSizeOptions.SEQ_LENGTH:
                    break

                image_chunk = images[:, i:i+TrainOptions.TRAIN_WINDOW]
                label_chunk = labels[:, i:i+TrainOptions.TRAIN_WINDOW]

                output_, hidden = transducer_model(image_chunk, hidden)

                loss = criterion(output_.contiguous().view(-1, num_classes), label_chunk.contiguous().view(-1))
                loss.backward()

                model_optimizer.step()
                total_loss += loss.item()
                total_images += image_chunk.size(0)
                hidden = hidden.detach()

            # update the progress bar
            avg_loss = (total_loss / total_images) if total_images else 0

            if train_mode is True and rank == 0:
                train_loss_logger.write(str(epoch + 1) + "," + str(batch_no) + "," + str(avg_loss) + "\n")

            if rank == 0:
                progress_bar.set_description("Loss: " + str(avg_loss))
                progress_bar.refresh()
                progress_bar.update(1)
                batch_no += 1

        if rank == 0:
            progress_bar.close()
        dist.barrier()

        if rank == 0:
            stats_dictioanry = test(test_file, batch_size, gpu_mode, transducer_model, num_workers,
                                    gru_layers, hidden_size, num_classes=ImageSizeOptions.TOTAL_LABELS)
            stats['loss'] = stats_dictioanry['loss']
            stats['accuracy'] = stats_dictioanry['accuracy']
            stats['loss_epoch'].append((epoch, stats_dictioanry['loss']))
            stats['accuracy_epoch'].append((epoch, stats_dictioanry['accuracy']))
        dist.barrier()

        # update the loggers
        if train_mode is True and rank == 0:
            # save the model after each epoch
            # encoder_model, decoder_model, encoder_optimizer, decoder_optimizer, hidden_size, layers, epoch,
            # file_name
            save_best_model(transducer_model, model_optimizer,
                            hidden_size, gru_layers, epoch, model_dir + "_epoch_" + str(epoch + 1) + '_checkpoint.pkl')

            test_loss_logger.write(str(epoch + 1) + "," + str(stats['loss']) + "," + str(stats['accuracy']) + "\n")
            confusion_matrix_logger.write(str(epoch + 1) + "\n" + str(stats_dictioanry['confusion_matrix']) + "\n")
            train_loss_logger.flush()
            test_loss_logger.flush()
            confusion_matrix_logger.flush()
        elif train_mode is False:
            # this setup is for hyperband
            if epoch + 1 >= 10 and stats['accuracy'] < 98:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: EARLY STOPPING AS THE MODEL NOT DOING WELL\n")
                return transducer_model, model_optimizer, stats

    if rank == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: Finished training\n")

    return transducer_model, model_optimizer, stats


def cleanup():
    dist.destroy_process_group()


def setup(rank, device_ids, args):
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = '12355'

    # initialize the process group
    dist.init_process_group("gloo", rank=rank, world_size=len(device_ids))

    train_file, test_file, batch_size, epochs, gpu_mode, num_workers, retrain_model, \
    retrain_model_path, gru_layers, hidden_size, learning_rate, weight_decay, model_dir, stats_dir, total_callers, \
    train_mode = args

    # issue with semaphore lock: https://github.com/pytorch/pytorch/issues/2517
    # mp.set_start_method('spawn')

    # Explicitly setting seed to make sure that models created in two processes
    # start from same random weights and biases. https://github.com/pytorch/pytorch/issues/2517
    torch.manual_seed(42)
    train(train_file, test_file, batch_size, epochs, gpu_mode, num_workers, retrain_model, retrain_model_path,
          gru_layers, hidden_size, learning_rate, weight_decay, model_dir, stats_dir, train_mode,
          total_callers, rank, device_ids[rank])
    cleanup()


def train_distributed(train_file, test_file, batch_size, epochs, gpu_mode, num_workers, retrain_model,
                          retrain_model_path, gru_layers, hidden_size, learning_rate, weight_decay, model_dir,
                          stats_dir, device_ids, total_callers, train_mode):

    args = (train_file, test_file, batch_size, epochs, gpu_mode, num_workers, retrain_model,
            retrain_model_path, gru_layers, hidden_size, learning_rate, weight_decay, model_dir,
            stats_dir, total_callers, train_mode)
    mp.spawn(setup,
             args=(device_ids, args),
             nprocs=len(device_ids),
             join=True)
