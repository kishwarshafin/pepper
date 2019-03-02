import sys
import torch
import torch.nn as nn
import os
from tqdm import tqdm
import numpy as np

# Custom generator for our dataset
from torch.utils.data import DataLoader
from modules.python.models.data_sampler import BalancedSampler
from modules.python.models.dataloader import SequenceDataset
from modules.python.TextColor import TextColor
from modules.python.models.ModelHander import ModelHandler
from modules.python.models.test import test
from modules.python.Options import ImageSizeOptions
"""
Train a model and return the model and optimizer trained.

Input:
- A train CSV containing training image set information (usually chr1-18)

Return:
- A trained model
"""
CLASS_WEIGHTS = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]


def save_best_model(encoder_model, decoder_model, encoder_optimizer, decoder_optimizer, hidden_size, layers, epoch,
                    file_name):
    """
    Save the best model
    :param encoder_model: A trained encoder model
    :param decoder_model: A trained decoder model
    :param encoder_optimizer: Encoder optimizer
    :param decoder_optimizer: Decoder optimizer
    :param file_name: Output file name
    :return:
    """
    if os.path.isfile(file_name):
        os.remove(file_name)
    ModelHandler.save_checkpoint({
        'encoder_state_dict': encoder_model.state_dict(),
        'decoder_state_dict': decoder_model.state_dict(),
        'encoder_optimizer': encoder_optimizer.state_dict(),
        'decoder_optimizer': decoder_optimizer.state_dict(),
        'hidden_size': hidden_size,
        'gru_layers': layers,
        'epochs': epoch,
    }, file_name)
    sys.stderr.write(TextColor.RED + "\nMODEL SAVED SUCCESSFULLY.\n" + TextColor.END)


def train(train_file, test_file, batch_size, epoch_limit, gpu_mode, num_workers, retrain_model,
          retrain_model_path, gru_layers, hidden_size, encoder_lr, encoder_decay, decoder_lr, decoder_decay,
          model_dir, stats_dir, train_mode):

    if train_mode is True:
        train_loss_logger = open(stats_dir + "train_loss.csv", 'w')
        test_loss_logger = open(stats_dir + "test_loss.csv", 'w')
        confusion_matrix_logger = open(stats_dir + "confusion_matrix.txt", 'w')
    else:
        train_loss_logger = None
        test_loss_logger = None
        confusion_matrix_logger = None

    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)
    train_data_set = SequenceDataset(train_file)
    train_loader = DataLoader(train_data_set,
                              batch_size=batch_size,
                              shuffle=True,
                              num_workers=num_workers,
                              pin_memory=gpu_mode)
    if retrain_model is True:
        if os.path.isfile(retrain_model_path) is False:
            sys.stderr.write(TextColor.RED + "ERROR: INVALID PATH TO RETRAIN PATH MODEL --retrain_model_path\n")
            exit(1)
        sys.stderr.write(TextColor.GREEN + "INFO: RETRAIN MODEL LOADING\n" + TextColor.END)
        encoder_model, decoder_model, hidden_size, gru_layers, prev_ite = \
            ModelHandler.load_model_for_training(retrain_model_path,
                                                 input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                 seq_len=ImageSizeOptions.SEQ_LENGTH,
                                                 num_classes=ImageSizeOptions.TOTAL_LABELS)

        if train_mode is True:
            epoch_limit = prev_ite + epoch_limit

        sys.stderr.write(TextColor.GREEN + "INFO: RETRAIN MODEL LOADED\n" + TextColor.END)
    else:
        encoder_model, decoder_model = ModelHandler.get_new_model(input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                                  gru_layers=gru_layers,
                                                                  hidden_size=hidden_size,
                                                                  seq_len=ImageSizeOptions.SEQ_LENGTH,
                                                                  num_classes=ImageSizeOptions.TOTAL_LABELS)
        prev_ite = 0

    encoder_optimizer = torch.optim.Adam(encoder_model.parameters(),
                                         lr=encoder_lr,
                                         weight_decay=encoder_decay)
    decoder_optimizer = torch.optim.Adam(decoder_model.parameters(),
                                         lr=decoder_lr,
                                         weight_decay=decoder_decay)

    if retrain_model is True:
        sys.stderr.write(TextColor.GREEN + "INFO: OPTIMIZER LOADING\n" + TextColor.END)
        encoder_optimizer, decoder_optimizer = ModelHandler.load_optimizer(encoder_optimizer, decoder_optimizer,
                                                                           retrain_model_path, gpu_mode)
        sys.stderr.write(TextColor.GREEN + "INFO: OPTIMIZER LOADED\n" + TextColor.END)

    if gpu_mode:
        encoder_model = torch.nn.DataParallel(encoder_model).cuda()
        decoder_model = torch.nn.DataParallel(decoder_model).cuda()

    class_weights = torch.FloatTensor(CLASS_WEIGHTS)
    # Loss
    criterion = nn.CrossEntropyLoss(weight=class_weights)

    if gpu_mode is True:
        criterion = criterion.cuda()

    start_epoch = prev_ite

    # Train the Model
    sys.stderr.write(TextColor.PURPLE + 'Training starting\n' + TextColor.END)
    stats = dict()
    stats['loss_epoch'] = []
    stats['accuracy_epoch'] = []
    sys.stderr.write(TextColor.PURPLE + 'Start: ' + str(start_epoch + 1) + ' End: ' + str(epoch_limit + 1) + "\n")
    for epoch in range(start_epoch, epoch_limit, 1):
        total_loss = 0
        total_images = 0
        sys.stderr.write(TextColor.BLUE + 'Train epoch: ' + str(epoch + 1) + "\n")
        # make sure the model is in train mode. BN is different in train and eval.
        encoder_model.train()
        decoder_model.train()
        batch_no = 1
        with tqdm(total=len(train_loader), desc='Loss', leave=True, ncols=100) as progress_bar:
            for images, labels in train_loader:
                # print(images.size(), labels.size())

                # from modules.python.helper.tensor_analyzer import analyze_tensor
                # for label in labels[0].data:
                #     print(label.item(), end='')
                # print()
                # analyze_tensor(images[0])
                # exit()
                if gpu_mode:
                    # encoder_hidden = encoder_hidden.cuda()
                    images = images.cuda()
                    labels = labels.cuda()

                encoder_hidden = torch.FloatTensor(images.size(0), gru_layers * 2, hidden_size).zero_()

                if gpu_mode:
                    encoder_hidden = encoder_hidden.cuda()

                encoder_optimizer.zero_grad()
                decoder_optimizer.zero_grad()

                loss = 0
                total_seq_length = images.size(2)
                start_index = ImageSizeOptions.CONTEXT_SIZE
                end_index = start_index + (ImageSizeOptions.SEQ_LENGTH - 2 * ImageSizeOptions.CONTEXT_SIZE)

                # from analysis.analyze_png_img import analyze_tensor
                # print(labels[0, :].data.numpy())
                # analyze_tensor(images[0, :, :, :])
                # exit()
                context_vector, hidden_encoder = encoder_model(images, encoder_hidden)

                for seq_index in range(start_index, end_index):
                    current_batch_size = images.size(0)
                    y = labels[:, seq_index - start_index]
                    attention_index = torch.from_numpy(np.asarray([seq_index] * current_batch_size)).view(-1, 1)

                    attention_index_onehot = torch.FloatTensor(current_batch_size, total_seq_length)

                    attention_index_onehot.zero_()
                    attention_index_onehot.scatter_(1, attention_index, 1)
                    # print("\n", seq_index, attention_index_onehot, y)
                    # exit()

                    output_dec, decoder_hidden, attn = decoder_model(attention_index_onehot,
                                                                     context_vector=context_vector,
                                                                     encoder_hidden=hidden_encoder)

                    # loss
                    loss += criterion(output_dec, y)

                loss.backward()
                encoder_optimizer.step()
                decoder_optimizer.step()

                total_loss += loss.item()
                total_images += labels.size(0)

                # update the progress bar
                avg_loss = (total_loss / total_images) if total_images else 0
                progress_bar.set_description("Loss: " + str(avg_loss))
                if train_mode is True:
                    train_loss_logger.write(str(epoch + 1) + "," + str(batch_no) + "," + str(avg_loss) + "\n")
                progress_bar.refresh()
                progress_bar.update(1)
                batch_no += 1

            progress_bar.close()

        stats_dictioanry = test(test_file, batch_size, gpu_mode, encoder_model, decoder_model, num_workers,
                                gru_layers, hidden_size, num_classes=ImageSizeOptions.TOTAL_LABELS)
        stats['loss'] = stats_dictioanry['loss']
        stats['accuracy'] = stats_dictioanry['accuracy']
        stats['loss_epoch'].append((epoch, stats_dictioanry['loss']))
        stats['accuracy_epoch'].append((epoch, stats_dictioanry['accuracy']))

        # update the loggers
        if train_mode is True:
            # save the model after each epoch
            # encoder_model, decoder_model, encoder_optimizer, decoder_optimizer, hidden_size, layers, epoch,
            # file_name
            save_best_model(encoder_model, decoder_model, encoder_optimizer, decoder_optimizer,
                            hidden_size, gru_layers, epoch, model_dir + "_epoch_" + str(epoch + 1) + '_checkpoint.pkl')

            test_loss_logger.write(str(epoch + 1) + "," + str(stats['loss']) + "," + str(stats['accuracy']) + "\n")
            confusion_matrix_logger.write(str(epoch + 1) + "\n" + str(stats_dictioanry['confusion_matrix']) + "\n")
            train_loss_logger.flush()
            test_loss_logger.flush()
            confusion_matrix_logger.flush()
        else:
            # this setup is for hyperband
            if epoch + 1 >= 2 and stats['accuracy'] < 90:
                sys.stderr.write(TextColor.PURPLE + 'EARLY STOPPING AS THE MODEL NOT DOING WELL\n' + TextColor.END)
                return encoder_model, decoder_model, encoder_optimizer, decoder_optimizer, stats

    sys.stderr.write(TextColor.PURPLE + 'Finished training\n' + TextColor.END)

    return encoder_model, decoder_model, encoder_optimizer, decoder_optimizer, stats

