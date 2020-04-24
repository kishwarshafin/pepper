import sys
import torch
import torch.nn as nn
import os
from tqdm import tqdm

# Custom generator for our dataset
from torch.utils.data import DataLoader
from pepper.modules.python.models.dataloader import SequenceDataset
from datetime import datetime
from pepper.modules.python.models.ModelHander import ModelHandler
from pepper.modules.python.models.test import test
from pepper.modules.python.Options import ImageSizeOptions, TrainOptions

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
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]  INFO: MODEL SAVED SUCCESSFULLY.\n")


def train(train_file, test_file, batch_size, epoch_limit, gpu_mode, num_workers, retrain_model,
          retrain_model_path, gru_layers, hidden_size, lr, decay, model_dir, stats_dir, train_mode):

    if train_mode is True:
        train_loss_logger = open(stats_dir + "train_loss.csv", 'w')
        test_loss_logger = open(stats_dir + "test_loss.csv", 'w')
        confusion_matrix_logger = open(stats_dir + "confusion_matrix.txt", 'w')
    else:
        train_loss_logger = None
        test_loss_logger = None
        confusion_matrix_logger = None

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: LOADING DATA\n")
    train_data_set = SequenceDataset(train_file)
    train_loader = DataLoader(train_data_set,
                              batch_size=batch_size,
                              shuffle=True,
                              num_workers=num_workers,
                              pin_memory=gpu_mode)
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
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL TRAINABLE PARAMETERS:\t" + str(param_count) + "\n")

    model_optimizer = torch.optim.Adam(transducer_model.parameters(), lr=lr, weight_decay=decay)
    lr_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(model_optimizer, 'min')

    if retrain_model is True:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: OPTIMIZER LOADING\n")
        model_optimizer = ModelHandler.load_simple_optimizer(model_optimizer, retrain_model_path, gpu_mode)
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: OPTIMIZER LOADED\n")

    if gpu_mode:
        transducer_model = transducer_model.cuda()

    class_weights = torch.Tensor(CLASS_WEIGHTS)
    # Loss
    criterion = nn.CrossEntropyLoss(class_weights)

    if gpu_mode is True:
        criterion = criterion.cuda()

    start_epoch = prev_ite

    # Train the Model
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: Training starting\n")
    stats = dict()
    stats['loss_epoch'] = []
    stats['accuracy_epoch'] = []
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] Start: " + str(start_epoch + 1) + " End: " + str(epoch_limit) + "\n")
    for epoch in range(start_epoch, epoch_limit, 1):
        total_loss = 0
        total_images = 0
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] Train epoch: " + str(epoch + 1) + "\n")
        # make sure the model is in train mode. BN is different in train and eval.

        batch_no = 1
        with tqdm(total=len(train_loader), desc='Loss', leave=True, ncols=100) as progress_bar:
            transducer_model.train()
            for images, labels in train_loader:
                labels = labels.type(torch.LongTensor)
                images = images.type(torch.FloatTensor)

                if gpu_mode:
                    # encoder_hidden = encoder_hidden.cuda()
                    images = images.cuda()
                    labels = labels.cuda()

                hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

                if gpu_mode:
                    hidden = hidden.cuda()

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
                progress_bar.set_description("Loss: " + str(avg_loss))

                if train_mode is True:
                    train_loss_logger.write(str(epoch + 1) + "," + str(batch_no) + "," + str(avg_loss) + "\n")
                progress_bar.refresh()
                progress_bar.update(1)
                batch_no += 1

            progress_bar.close()

        stats_dictioanry = test(test_file, batch_size, gpu_mode, transducer_model, num_workers,
                                gru_layers, hidden_size, num_classes=ImageSizeOptions.TOTAL_LABELS)
        stats['loss'] = stats_dictioanry['loss']
        stats['accuracy'] = stats_dictioanry['accuracy']
        stats['loss_epoch'].append((epoch, stats_dictioanry['loss']))
        stats['accuracy_epoch'].append((epoch, stats_dictioanry['accuracy']))

        lr_scheduler.step(stats['loss'])

        # update the loggers
        if train_mode is True:
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
        else:
            # this setup is for hyperband
            if epoch + 1 >= 10 and stats['accuracy'] < 98:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: EARLY STOPPING AS THE MODEL NOT DOING WELL\n")
                return transducer_model, model_optimizer, stats

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: Finished training\n")

    return transducer_model, model_optimizer, stats

