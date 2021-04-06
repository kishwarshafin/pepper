import sys
import torch
import torchnet.meter as meter
import torch.nn as nn
import time
from datetime import datetime
from torch.utils.data import DataLoader
from pepper_variant.modules.python.models.dataloader import SequenceDatasetHP
from pepper_variant.modules.python.Options import ImageSizeOptionsHP, TrainOptions
"""
This script will evaluate a model and return the loss value.

Input:
- A trained model
- A test CSV file to evaluate

Returns:
- Loss value
"""
label_decoder = {0: '*', 1: 'A', 2: 'C', 3: 'G', 4: 'T', 5: '#'}


def test_hp(data_file, batch_size, gpu_mode, transducer_model, num_workers, gru_layers, hidden_size,
            num_classes=ImageSizeOptionsHP.TOTAL_LABELS, print_details=False):
    start_time = time.time()
    # data loader
    test_data = SequenceDatasetHP(data_file)
    test_loader = DataLoader(test_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers,
                             pin_memory=gpu_mode)

    # set the evaluation mode of the model
    transducer_model.eval()

    # Loss
    criterion = nn.CrossEntropyLoss()

    if gpu_mode is True:
        criterion = criterion.cuda()

    # Test the Model
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST STARTING\n")
    confusion_matrix = meter.ConfusionMeter(num_classes)

    total_loss = 0
    total_images = 0
    accuracy = 0

    with torch.no_grad():
        for ii, (images, labels) in enumerate(test_loader):
            labels = labels.type(torch.LongTensor)
            images = images.type(torch.FloatTensor)
            if gpu_mode:
                # encoder_hidden = encoder_hidden.cuda()
                images = images.cuda()
                labels = labels.cuda()

            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            if gpu_mode:
                hidden = hidden.cuda()

            for i in range(0, ImageSizeOptionsHP.SEQ_LENGTH, TrainOptions.WINDOW_JUMP):
                if i + TrainOptions.TRAIN_WINDOW > ImageSizeOptionsHP.SEQ_LENGTH:
                    break

                image_chunk = images[:, i:i+TrainOptions.TRAIN_WINDOW]
                label_chunk = labels[:, i:i+TrainOptions.TRAIN_WINDOW]
                output_, hidden = transducer_model(image_chunk, hidden)

                loss = criterion(output_.contiguous().view(-1, num_classes), label_chunk.contiguous().view(-1))

                confusion_matrix.add(output_.data.contiguous().view(-1, num_classes),
                                     label_chunk.data.contiguous().view(-1))

                total_loss += loss.item()
                total_images += images.size(0)

            if (ii + 1) % 10 == 0:
                cm_value = confusion_matrix.value()
                denom = cm_value.sum() if cm_value.sum() > 0 else 1.0
                accuracy = 100.0 * (cm_value[0][0] + cm_value[1][1] + cm_value[2][2]
                                    + cm_value[3][3] + cm_value[4][4]) / denom
                percent_complete = int((100 * (ii+1)) / len(test_loader))
                time_now = time.time()
                mins = int((time_now - start_time) / 60)
                secs = int((time_now - start_time)) % 60

                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO:"
                                 + " BATCH: " + str(ii+1) + "/" + str(len(test_loader))
                                 + " ACCURACY: " + str(round(accuracy, 5))
                                 + " COMPLETE (" + str(percent_complete) + "%)"
                                 + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                sys.stderr.flush()

    avg_loss = total_loss / total_images if total_images else 0

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST LOSS: " + str(avg_loss) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CONFUSION MATRIX: \n" + str(confusion_matrix.conf) + "\n")
    sys.stderr.flush()

    return {'loss': avg_loss, 'accuracy': accuracy, 'confusion_matrix': str(confusion_matrix.conf)}
