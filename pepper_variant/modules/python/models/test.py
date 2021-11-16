import sys
import torch
import torchnet.meter as meter
import torch.nn as nn
import time
from datetime import datetime
from torch.utils.data import DataLoader
from pepper_variant.modules.python.models.dataloader import SequenceDataset
from pepper_variant.modules.python.Options import ImageSizeOptions, TrainOptions
"""
This script will evaluate a model and return the loss value.

Input:
- A trained model
- A test CSV file to evaluate

Returns:
- Loss value
"""


def test(data_file, batch_size, gpu_mode, transducer_model, num_workers, gru_layers, hidden_size, train_with_base, num_classes=ImageSizeOptions.TOTAL_LABELS):
    start_time = time.time()
    # data loader
    test_data = SequenceDataset(data_file)
    test_loader = DataLoader(test_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=0,
                             pin_memory=False,
                             drop_last=True)

    # set the evaluation mode of the model
    transducer_model.eval()

    # Loss
    if train_with_base:
        criterion_base = nn.CrossEntropyLoss()
    criterion_type = nn.CrossEntropyLoss()

    if gpu_mode is True:
        if train_with_base:
            criterion_base = criterion_base.cuda()
        criterion_type = criterion_type.cuda()

    # Test the Model
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST STARTING\n")
    sys.stderr.flush()
    confusion_matrix_base = meter.ConfusionMeter(ImageSizeOptions.TOTAL_LABELS)
    confusion_matrix_type = meter.ConfusionMeter(ImageSizeOptions.TOTAL_TYPE_LABELS)

    total_loss = 0
    total_base_loss = 0
    total_type_loss = 0
    total_images = 0

    with torch.no_grad():
        for ii, (images, base_labels, type_labels) in enumerate(test_loader):
            base_labels = base_labels.type(torch.LongTensor)
            type_labels = type_labels.type(torch.LongTensor)
            images = images.type(torch.FloatTensor)
            if gpu_mode:
                images = images.cuda()
                base_labels = base_labels.cuda()
                type_labels = type_labels.cuda()


            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
            cell_state = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            if gpu_mode:
                cell_state = cell_state.cuda()
                hidden = hidden.cuda()

            output_base, output_type = transducer_model(images, hidden, cell_state, train_mode=True)

            if train_with_base:
                loss_base = criterion_base(output_base.contiguous().view(-1, ImageSizeOptions.TOTAL_LABELS), base_labels.contiguous().view(-1))
            loss_type = criterion_type(output_type.contiguous().view(-1, ImageSizeOptions.TOTAL_TYPE_LABELS), type_labels.contiguous().view(-1))

            if train_with_base:
                loss = 0.3 * loss_base + loss_type
            else:
                loss = loss_type

            if train_with_base:
                confusion_matrix_base.add(output_base.data.contiguous().view(-1, num_classes), base_labels.data.contiguous().view(-1))

            confusion_matrix_type.add(output_type.data.contiguous().view(-1, ImageSizeOptions.TOTAL_TYPE_LABELS), type_labels.data.contiguous().view(-1))

            total_loss += loss.item()
            if train_with_base:
                total_base_loss += loss_base.item()
            total_type_loss += loss_type.item()
            total_images += images.size(0)

            if (ii + 1) % 100 == 0:
                if train_with_base:
                    cm_value_base = confusion_matrix_base.value()
                    denom_base = cm_value_base.sum() if cm_value_base.sum() > 0 else 1.0

                    total_accurate_base = 0
                    for i in range(0, ImageSizeOptions.TOTAL_LABELS):
                        total_accurate_base = total_accurate_base + cm_value_base[i][i]

                    accuracy_base = (100.0 * total_accurate_base) / denom_base

                cm_value_type = confusion_matrix_type.value()
                denom_type = cm_value_type.sum() if cm_value_type.sum() > 0 else 1.0

                total_accurate_type = 0
                for i in range(0, ImageSizeOptions.TOTAL_TYPE_LABELS):
                    total_accurate_type = total_accurate_type + cm_value_type[i][i]

                accuracy_type = (100.0 * total_accurate_type) / denom_type

                percent_complete = int((100 * (ii+1)) / len(test_loader))
                time_now = time.time()
                mins = int((time_now - start_time) / 60)
                secs = int((time_now - start_time)) % 60

                if train_with_base:
                    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO:"
                                     + " BATCH: " + str(ii+1) + "/" + str(len(test_loader))
                                     + " ACCURACY BASE: " + str(round(accuracy_base, 5))
                                     + " ACCURACY TYPE: " + str(round(accuracy_type, 5))
                                     + " COMPLETE (" + str(percent_complete) + "%)"
                                     + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                else:
                    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO:"
                                     + " BATCH: " + str(ii+1) + "/" + str(len(test_loader))
                                     + " ACCURACY TYPE: " + str(round(accuracy_type, 5))
                                     + " COMPLETE (" + str(percent_complete) + "%)"
                                     + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                sys.stderr.flush()
    avg_loss = total_loss / total_images if total_images else 0
    avg_base_loss = total_base_loss / total_images if total_images else 0
    avg_type_loss = total_type_loss / total_images if total_images else 0

    # base confusion matrix
    if train_with_base:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST BASE LOSS: " + str(avg_base_loss) + "\n")
        sys.stderr.write("Confusion Matrix BASE:" + "\n")
        sys.stderr.write("         ")
        sys.stderr.flush()
        for label in ImageSizeOptions.decoded_base_labels:
            sys.stderr.write(str(label) + '      ')
        sys.stderr.write("\n")
        sys.stderr.flush()

        for i, row in enumerate(confusion_matrix_base.value()):
            sys.stderr.write(str(ImageSizeOptions.decoded_base_labels[i]) + '   ')
            for j, val in enumerate(row):
                sys.stderr.write("{0:6d}".format(val) + '  ')
            sys.stderr.write("\n")
        sys.stderr.flush()

    # type Confusion matrix
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST TYPE LOSS: " + str(avg_type_loss) + "\n")
    sys.stderr.write("Confusion Matrix TYPE:" + "\n")
    sys.stderr.write("            ")
    sys.stderr.flush()
    for label in ImageSizeOptions.decoded_labels:
        sys.stderr.write(str(label) + '    ')
    sys.stderr.write("\n")
    sys.stderr.flush()

    for i, row in enumerate(confusion_matrix_type.value()):
        sys.stderr.write(str(ImageSizeOptions.decoded_labels[i]) + '   ')
        for j, val in enumerate(row):
            sys.stderr.write("{0:9d}".format(val) + '  ')
        sys.stderr.write("\n")
    sys.stderr.flush()

    cm_value = confusion_matrix_base.value()
    denom = cm_value.sum() if cm_value.sum() > 0 else 1.0

    total_type_accurate = 0
    for i in range(0, ImageSizeOptions.TOTAL_TYPE_LABELS):
        total_type_accurate = total_type_accurate + cm_value[i][i]

    base_accuracy = (100.0 * total_type_accurate) / denom

    return {'loss': avg_loss, 'base_accuracy': base_accuracy, 'base_confusion_matrix': confusion_matrix_base, 'type_confusion_matrix': confusion_matrix_type}
