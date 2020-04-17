import sys
import torch
import time
import torchnet.meter as meter
import torch.nn as nn
from torch.utils.data import DataLoader
from datetime import datetime
from pepper_snp.modules.python.models.dataloader import SequenceDataset
from pepper_snp.modules.python.Options import ImageSizeOptions, TrainOptions
"""
This script will evaluate a model and return the loss value.

Input:
- A trained model
- A test CSV file to evaluate

Returns:
- Loss value
"""
CLASS_WEIGHTS = [0.3, 1.0, 1.0, 1.0, 1.0]
label_decoder = {0: '*', 1: 'A', 2: 'C', 3: 'G', 4: 'T', 5: '#'}


def test(data_file, batch_size, gpu_mode, transducer_model, num_workers, gru_layers, hidden_size,
         num_classes=ImageSizeOptions.TOTAL_LABELS):
    start_time = time.time()

    # data loader
    test_data = SequenceDataset(data_file)
    test_loader = DataLoader(test_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers,
                             pin_memory=gpu_mode)
    # sys.stderr.write(TextColor.CYAN + 'Test data loaded\n')

    # set the evaluation mode of the model
    transducer_model.eval()
    class_weights = torch.Tensor(ImageSizeOptions.class_weights)
    # Loss
    criterion = nn.CrossEntropyLoss(class_weights)

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
                images = images.cuda()
                labels = labels.cuda()

            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            if gpu_mode:
                hidden = hidden.cuda()

            for i in range(0, ImageSizeOptions.SEQ_LENGTH, TrainOptions.WINDOW_JUMP):
                if i + TrainOptions.TRAIN_WINDOW > ImageSizeOptions.SEQ_LENGTH:
                    break

                image_chunk = images[:, i:i+TrainOptions.TRAIN_WINDOW]
                label_chunk = labels[:, i:i+TrainOptions.TRAIN_WINDOW]
                output_, hidden = transducer_model(image_chunk, hidden)

                loss = criterion(output_.contiguous().view(-1, num_classes), label_chunk.contiguous().view(-1))

                confusion_matrix.add(output_.data.contiguous().view(-1, num_classes),
                                     label_chunk.data.contiguous().view(-1))

                total_loss += loss.item()
                total_images += images.size(0)

            cm_value = confusion_matrix.value()
            denom = cm_value.sum() if cm_value.sum() > 0 else 1.0

            total_accurate = 0
            for i in range(0, ImageSizeOptions.TOTAL_LABELS):
                total_accurate = total_accurate + cm_value[i][i]

            accuracy = (100.0 * total_accurate) / denom

            if (ii + 1) % 10 == 0:
                percent_complete = int((100 * (ii+1)) / len(test_loader))
                time_now = time.time()
                mins = int((time_now - start_time) / 60)
                secs = int((time_now - start_time)) % 60

                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO:"
                                 + " BATCH: " + str(ii+1) + "/" + str(len(test_loader))
                                 + " ACCURACY: " + str(round(accuracy, 5))
                                 + " COMPLETE (" + str(percent_complete) + "%)"
                                 + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")

    avg_loss = total_loss / total_images if total_images else 0

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST LOSS: " + str(avg_loss) + "\n")

    sys.stderr.write("Confusion Matrix:" + "\n")
    sys.stderr.write("            ")
    for label in ImageSizeOptions.decoded_labels:
        sys.stderr.write(str(label) + '         ')
    sys.stderr.write("\n")

    for i, row in enumerate(confusion_matrix.value()):
        sys.stderr.write(str(ImageSizeOptions.decoded_labels[i]) + '   ')
        for j, val in enumerate(row):
            sys.stderr.write("{0:9d}".format(val) + '  ')
        sys.stderr.write("\n")

    return {'loss': avg_loss, 'accuracy': accuracy, 'confusion_matrix': str(confusion_matrix.conf)}
