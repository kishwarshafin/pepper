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


def test(data_file, batch_size, gpu_mode, transducer_model, num_workers, gru_layers, hidden_size,
         num_classes=ImageSizeOptions.TOTAL_LABELS, num_type_classes=ImageSizeOptions.TOTAL_TYPE_LABELS, print_details=False):
    start_time = time.time()
    # data loader
    test_data = SequenceDataset(data_file)
    test_loader = DataLoader(test_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers,
                             pin_memory=gpu_mode)

    # set the evaluation mode of the model
    transducer_model.eval()
    class_weights = torch.Tensor(ImageSizeOptions.class_weights)
    class_weights_type = torch.Tensor(ImageSizeOptions.class_weights_type)
    # Loss
    criterion_base = nn.NLLLoss(class_weights)
    criterion_type = nn.NLLLoss(class_weights_type)

    if gpu_mode is True:
        criterion_base = criterion_base.cuda()
        criterion_type = criterion_type.cuda()

    # Test the Model
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST STARTING\n")
    confusion_matrix = meter.ConfusionMeter(num_classes)
    confusion_matrix_type = meter.ConfusionMeter(num_type_classes)

    total_loss = 0
    total_base_loss = 0
    total_type_loss = 0
    total_images = 0
    accuracy = 0

    with torch.no_grad():
        for ii, (images, labels, type_labels) in enumerate(test_loader):
            labels = labels.type(torch.LongTensor)
            type_labels = type_labels.type(torch.LongTensor)
            images = images.type(torch.FloatTensor)
            if gpu_mode:
                # encoder_hidden = encoder_hidden.cuda()
                images = images.cuda()
                labels = labels.cuda()
                type_labels = type_labels.cuda()

            hidden = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)
            cell_state = torch.zeros(images.size(0), 2 * TrainOptions.GRU_LAYERS, TrainOptions.HIDDEN_SIZE)

            if gpu_mode:
                hidden = hidden.cuda()

            output_base, output_type = transducer_model(images, hidden, cell_state, train_mode=True)
            # output_base = transducer_model(images, hidden, cell_state, train_mode=True)

            loss_base = criterion_base(output_base.contiguous().view(-1, num_classes), labels.contiguous().view(-1))
            loss_type = criterion_type(output_type.contiguous().view(-1, num_type_classes), type_labels.contiguous().view(-1))

            loss = loss_base + loss_type
            confusion_matrix.add(output_base.data.contiguous().view(-1, num_classes),
                                 labels.data.contiguous().view(-1))

            confusion_matrix_type.add(output_type.data.contiguous().view(-1, num_type_classes),
                                      type_labels.data.contiguous().view(-1))

            total_loss += loss.item()
            total_base_loss += loss_base.item()
            total_type_loss += loss_type.item()
            total_images += images.size(0)

            if (ii + 1) % 50 == 0:
                cm_value = confusion_matrix.value()
                denom = cm_value.sum() if cm_value.sum() > 0 else 1.0

                total_accurate = 0
                for i in range(0, ImageSizeOptions.TOTAL_LABELS):
                    total_accurate = total_accurate + cm_value[i][i]

                accuracy = (100.0 * total_accurate) / denom

                cm_value_type = confusion_matrix_type.value()
                denom_type = cm_value_type.sum() if cm_value_type.sum() > 0 else 1.0

                for i in range(0, ImageSizeOptions.TOTAL_TYPE_LABELS):
                    denom_type = denom_type - cm_value_type[0][i]

                total_accurate_type = 0
                for i in range(1, ImageSizeOptions.TOTAL_TYPE_LABELS):
                    total_accurate_type = total_accurate_type + cm_value_type[i][i]
                accuracy_type = (100.0 * total_accurate_type) / denom_type

                percent_complete = int((100 * (ii+1)) / len(test_loader))
                time_now = time.time()
                mins = int((time_now - start_time) / 60)
                secs = int((time_now - start_time)) % 60

                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO:"
                                 + " BATCH: " + str(ii+1) + "/" + str(len(test_loader))
                                 + " ACCURACY BASE: " + str(round(accuracy, 5))
                                 + " ACCURACY TYPE: " + str(round(accuracy_type, 5))
                                 + " COMPLETE (" + str(percent_complete) + "%)"
                                 + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                sys.stderr.flush()

    avg_loss = total_loss / total_images if total_images else 0
    avg_base_loss = total_base_loss / total_images if total_images else 0
    avg_type_loss = total_type_loss / total_images if total_images else 0

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
    sys.stderr.flush()

    sys.stderr.write("Type Confusion Matrix:" + "\n")
    sys.stderr.write("            ")
    for label in ImageSizeOptions.decoded_type_labels:
        sys.stderr.write(str(label) + '         ')
    sys.stderr.write("\n")

    for i, row in enumerate(confusion_matrix_type.value()):
        sys.stderr.write(str(ImageSizeOptions.decoded_type_labels[i]) + '   ')
        for j, val in enumerate(row):
            sys.stderr.write("{0:9d}".format(val) + '  ')
        sys.stderr.write("\n")
    sys.stderr.flush()

    cm_value = confusion_matrix.value()
    denom = cm_value.sum() if cm_value.sum() > 0 else 1.0

    total_base_accurate = 0
    for i in range(0, ImageSizeOptions.TOTAL_LABELS):
        total_base_accurate = total_base_accurate + cm_value[i][i]

    base_accuracy = (100.0 * total_base_accurate) / denom

    cm_value_type = confusion_matrix_type.value()
    denom_type = cm_value_type.sum() if cm_value_type.sum() > 0 else 1.0

    for i in range(0, ImageSizeOptions.TOTAL_TYPE_LABELS):
        denom_type = denom_type - cm_value_type[0][i]

    total_accurate_type = 0
    for i in range(1, ImageSizeOptions.TOTAL_TYPE_LABELS):
        total_accurate_type = total_accurate_type + cm_value_type[i][i]
    type_accuracy = (100.0 * total_accurate_type) / denom_type

    return {'loss': avg_loss, 'base_loss': avg_base_loss, 'type_loss': avg_type_loss, 'base_accuracy': base_accuracy, 'type_accuracy': type_accuracy, 'base_confusion_matrix': confusion_matrix, 'type_confusion_matrix': confusion_matrix_type}
