import sys
import torch
import os
import h5py
from os.path import isfile, join
from os import listdir
import time
import torch.distributed as dist
import torch.nn as nn
import torch.multiprocessing as mp
from datetime import datetime
# Custom generator for our dataset
from torch.utils.data import DataLoader
from pepper_variant.modules.python.models.dataloader import SequenceDataset
from pepper_variant.modules.python.models.ModelHander import ModelHandler
from pepper_variant.modules.python.models.test import test
from pepper_variant.modules.python.Options import ImageSizeOptions, TrainOptions

os.environ['PYTHONWARNINGS'] = 'ignore:semaphore_tracker:UserWarning'

"""
Train a model and return the model and optimizer trained.

Input:
- A train CSV containing training image set information (usually chr1-18)

Return:
- A trained model
"""


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
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL" + file_name + " SAVED SUCCESSFULLY.\n")
    sys.stderr.flush()


def train(train_file, test_file, batch_size, test_batch_size, step_size, epoch_limit, gpu_mode, num_workers, retrain_model,
          retrain_model_path, gru_layers, hidden_size, lr, decay, model_dir, stats_dir, train_mode,
          world_size, rank, device_id):

    if train_mode is True and rank == 0:
        train_loss_logger = open(stats_dir + "train_loss.csv", 'w')
        test_loss_logger = open(stats_dir + "test_loss.csv", 'w')
        confusion_matrix_logger = open(stats_dir + "base_confusion_matrix.txt", 'w')
    else:
        train_loss_logger = None
        test_loss_logger = None
        confusion_matrix_logger = None

    torch.cuda.set_device(device_id)

    if rank == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: LOADING DATA\n")
        sys.stderr.flush()

    train_data_set = SequenceDataset(train_file)

    train_sampler = torch.utils.data.distributed.DistributedSampler(
        train_data_set,
        num_replicas=world_size,
        rank=rank
    )

    # shuffle is off because we are using distributed sampler, which has shuffle on
    train_loader = torch.utils.data.DataLoader(
        dataset=train_data_set,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0,
        pin_memory=True,
        sampler=train_sampler)

    num_classes = ImageSizeOptions.TOTAL_LABELS
    num_type_classes = ImageSizeOptions.TOTAL_TYPE_LABELS

    if retrain_model is True:
        if os.path.isfile(retrain_model_path) is False:
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: INVALID PATH TO RETRAIN PATH MODEL --retrain_model_path\n")
            exit(1)
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RETRAIN MODEL LOADING\n")
        transducer_model, hidden_size, gru_layers, prev_ite = \
            ModelHandler.load_simple_model_for_training(retrain_model,
                                                        image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                        num_classes=ImageSizeOptions.TOTAL_LABELS,
                                                        num_type_classes=ImageSizeOptions.TOTAL_TYPE_LABELS)

        if train_mode is True:
            epoch_limit = prev_ite + epoch_limit

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RETRAIN MODEL LOADED\n")
    else:
        transducer_model = ModelHandler.get_new_gru_model(image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                          gru_layers=gru_layers,
                                                          hidden_size=hidden_size,
                                                          num_classes=num_classes,
                                                          num_classes_type=num_type_classes)
        prev_ite = 0

    param_count = sum(p.numel() for p in transducer_model.parameters() if p.requires_grad)
    if rank == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL TRAINABLE PARAMETERS:\t" + str(param_count) + "\n")
        sys.stderr.flush()

    model_optimizer = torch.optim.Adam(transducer_model.parameters(), lr=lr, weight_decay=decay)

    if retrain_model is True:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: OPTIMIZER LOADING\n")
        model_optimizer = ModelHandler.load_simple_optimizer(model_optimizer, retrain_model_path, gpu_mode)
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: OPTIMIZER LOADED\n")
        sys.stderr.flush()

    if gpu_mode:
        transducer_model = transducer_model.to(device_id)
        transducer_model = nn.parallel.DistributedDataParallel(transducer_model, device_ids=[device_id])

    class_weights = torch.tensor(ImageSizeOptions.class_weights).to(device_id)
    # class_weights_type = torch.tensor(ImageSizeOptions.class_weights_type).to(device_id)
    # Loss
    criterion_base = nn.NLLLoss(weight=class_weights)
    # criterion_type = nn.NLLLoss(weight=class_weights_type)

    if gpu_mode is True:
        criterion_base = criterion_base.to(device_id)
        # criterion_type = criterion_type.to(device_id)

    start_epoch = prev_ite

    # Train the Model
    if rank == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TRAINING STARTING\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: START: " + str(start_epoch + 1) + " END: " + str(epoch_limit) + "\n")
        sys.stderr.flush()

    # Creates a GradScaler once at the beginning of training.
    # scaler = torch.cuda.amp.GradScaler()

    start_iteration = 0
    if step_size == 0:
        step_size = len(train_loader)

    iteration_limit = int((step_size * epoch_limit) / len(train_loader)) + 1
    step_no = 0
    epoch = 0

    if rank == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ITERATION LIMIT: " + str(iteration_limit) + "\n")
        sys.stderr.flush()

    for iteration in range(start_iteration, iteration_limit, 1):
        train_sampler.set_epoch(iteration)
        start_time = time.time()
        total_loss = 0
        total_base_loss = 0
        total_type_loss = 0
        total_images = 0

        for images, labels, type_labels in train_loader:
            # make sure the model is in train mode.
            transducer_model = transducer_model.train()

            # image and label handling
            labels = labels.type(torch.LongTensor)
            # type_labels = type_labels.type(torch.LongTensor)
            images = images.type(torch.FloatTensor)
            if gpu_mode:
                images = images.to(device_id)
                labels = labels.to(device_id)
                # type_labels = type_labels.to(device_id)

            # generate hidden and cell state
            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
            cell_state = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            if gpu_mode:
                hidden = hidden.to(device_id)
                cell_state = cell_state.to(device_id)

            # set the optimizer to zero grad
            model_optimizer.zero_grad()

            # Runs the forward pass with autocasting.
            # with torch.cuda.amp.autocast():
            # turned of mixed precision training

            output_base = transducer_model(images, hidden, cell_state, train_mode)
            # output_base, output_type = transducer_model(images, hidden, cell_state, train_mode)

            loss = criterion_base(output_base.contiguous().view(-1, num_classes), labels.contiguous().view(-1))
            # loss_type = criterion_type(output_type.contiguous().view(-1, num_type_classes), type_labels.contiguous().view(-1))
            # loss = loss_base + loss_type

            # Scales loss.
            # scaler.scale(loss).backward()
            loss.backward()

            # scaler.step() first unscales the gradients of the optimizer's assigned params.
            # scaler.step(model_optimizer)
            model_optimizer.step()

            # Updates the scale for next iteration.
            # scaler.update()

            # done calculating the loss
            total_loss += loss.item()
            # total_type_loss += loss_type.item()
            total_loss += loss.item()
            total_images += images.size(0)

            # update the progress bar
            avg_loss = (total_loss / total_images) if total_images else 0
            # avg_base_loss = (total_base_loss / total_images) if total_images else 0
            # avg_type_loss = (total_type_loss / total_images) if total_images else 0

            if rank == 0 and step_no % 10 == 0 and step_no > 0:
                percent_complete = int((float(step_no) / float((iteration + 1) * len(train_loader))) * 100)
                time_now = time.time()
                mins = int((time_now - start_time) / 60)
                secs = int((time_now - start_time)) % 60

                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: "
                                 + " ITERATION: " + str(iteration + 1)
                                 + " EPOCH: " + str(epoch + 1)
                                 + " STEP: " + str(step_no) + "/" + str((epoch + 1) * step_size) + "/" + str((iteration + 1) * len(train_loader))
                                 + " LOSS: " + "{:.9f}".format(avg_loss)
                                 # + " BASE LOSS: " + "{:.9f}".format(avg_base_loss)
                                 # + " TYPE LOSS: " + "{:.9f}".format(avg_type_loss)
                                 + " COMPLETE (" + str(percent_complete) + "%)"
                                 + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                sys.stderr.flush()

            if step_no % step_size == 0 and step_no > 0:
                dist.barrier()

                if rank == 0:
                    transducer_model = transducer_model.eval()
                    torch.cuda.empty_cache()

                    stats_dictioanry = test(test_file, test_batch_size, gpu_mode, transducer_model, num_workers,
                                            gru_layers, hidden_size, num_classes=ImageSizeOptions.TOTAL_LABELS, num_type_classes=ImageSizeOptions.TOTAL_TYPE_LABELS)

                    save_best_model(transducer_model, model_optimizer, hidden_size, gru_layers, epoch, model_dir + "PEPPER_VARIANT_STEP_" + str(step_no) + '_checkpoint.pkl')

                    train_loss_logger.write(str(step_no) + "," + str(avg_loss) + "\n")
                    # test_loss_logger.write(str(step_no) + "," + str(stats_dictioanry['loss']) + "," + str(stats_dictioanry['base_loss']) + "," + str(stats_dictioanry['type_loss']) + "," + str(stats_dictioanry['base_accuracy']) + "," + str(stats_dictioanry['type_accuracy']) + "\n")
                    test_loss_logger.write(str(step_no) + "," + str(stats_dictioanry['loss']) + "," + str(stats_dictioanry['base_accuracy']) + "\n")
                    confusion_matrix_logger.write(str(step_no) + "\n")

                    # print confusion matrix to file
                    confusion_matrix_logger.write("Confusion Matrix:" + "\n")
                    confusion_matrix_logger.write("            ")
                    for label in ImageSizeOptions.decoded_labels:
                        confusion_matrix_logger.write(str(label) + '         ')
                    confusion_matrix_logger.write("\n")

                    for i, row in enumerate(stats_dictioanry['base_confusion_matrix'].value()):
                        confusion_matrix_logger.write(str(ImageSizeOptions.decoded_labels[i]) + '   ')
                        for j, val in enumerate(row):
                            confusion_matrix_logger.write("{0:9d}".format(val) + '  ')
                        confusion_matrix_logger.write("\n")
                    confusion_matrix_logger.flush()

                    # confusion_matrix_logger.write("Type Confusion Matrix:" + "\n")
                    # confusion_matrix_logger.write("            ")
                    # for label in ImageSizeOptions.decoded_type_labels:
                    #     confusion_matrix_logger.write(str(label) + '         ')
                    # confusion_matrix_logger.write("\n")
                    #
                    # for i, row in enumerate(stats_dictioanry['type_confusion_matrix'].value()):
                    #     confusion_matrix_logger.write(str(ImageSizeOptions.decoded_type_labels[i]) + '   ')
                    #     for j, val in enumerate(row):
                    #         confusion_matrix_logger.write("{0:9d}".format(val) + '  ')
                    #     confusion_matrix_logger.write("\n")
                    # confusion_matrix_logger.flush()

                    train_loss_logger.flush()
                    test_loss_logger.flush()
                    confusion_matrix_logger.flush()

                    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST COMPLETED.\n")
                    sys.stderr.flush()
                    transducer_model = transducer_model.train()

                epoch += 1

                dist.barrier()

            # increase step size
            step_no += 1

            if epoch == epoch_limit:
                dist.barrier()
                break

        if rank == 0:
            time_now = time.time()
            mins = int((time_now - start_time) / 60)
            secs = int((time_now - start_time)) % 60
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ELAPSED TIME FOR ONE ITERATION: " + str(mins) + " Min " + str(secs) + " Sec\n")

        if epoch == epoch_limit:
            dist.barrier()
            break

    if rank == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED TRAINING\n")


def cleanup():
    dist.destroy_process_group()


def setup(rank, device_ids, args):
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = '12355'

    # initialize the process group
    dist.init_process_group("gloo", rank=rank, world_size=len(device_ids))

    train_file, test_file, batch_size, test_batch_size, step_size, epochs, gpu_mode, num_workers, retrain_model, \
    retrain_model_path, gru_layers, hidden_size, learning_rate, weight_decay, model_dir, stats_dir, total_callers, \
    train_mode = args

    # issue with semaphore lock: https://github.com/pytorch/pytorch/issues/2517
    # mp.set_start_method('spawn')

    # Explicitly setting seed to make sure that models created in two processes
    # start from same random weights and biases. https://github.com/pytorch/pytorch/issues/2517
    torch.manual_seed(42)
    train(train_file, test_file, batch_size, test_batch_size, step_size, epochs, gpu_mode, num_workers, retrain_model, retrain_model_path,
          gru_layers, hidden_size, learning_rate, weight_decay, model_dir, stats_dir, train_mode,
          total_callers, rank, device_ids[rank])
    cleanup()

def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-3:] == 'hdf']
    return file_paths


def train_distributed(train_file, test_file, batch_size, test_batch_size, step_size, epochs, gpu_mode, num_workers, retrain_model,
                      retrain_model_path, gru_layers, hidden_size, learning_rate, weight_decay, model_dir,
                      stats_dir, device_ids, total_callers, train_mode):

    args = (train_file, test_file, batch_size, test_batch_size, step_size, epochs, gpu_mode, num_workers, retrain_model,
            retrain_model_path, gru_layers, hidden_size, learning_rate, weight_decay, model_dir,
            stats_dir, total_callers, train_mode)

    mp.spawn(setup,
             args=(device_ids, args),
             nprocs=len(device_ids),
             join=True)
