import sys
import torch
import torchnet.meter as meter
import torch.nn as nn
import time
from datetime import datetime
from torch.utils.data import DataLoader
from pepper_variant.modules.python.models.dataloader import SequenceDataset
from pepper_variant.modules.python.Options import ImageSizeOptions, TrainOptions


def test(data_file, batch_size, gpu_mode, transducer_model, num_workers):
    start_time = time.time()
    # data loader
    test_data = SequenceDataset(data_file)
    test_loader = DataLoader(test_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers,
                             pin_memory=False,
                             drop_last=True)

    # set the evaluation mode of the model
    transducer_model.eval()

    # Loss
    criterion_type = nn.CrossEntropyLoss(reduction='sum')

    if gpu_mode is True:
        criterion_type = criterion_type.cuda()

    # Test the Model
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST STARTING\n")
    sys.stderr.flush()
    confusion_matrix = meter.ConfusionMeter(ImageSizeOptions.TOTAL_TYPE_LABELS)

    total_loss = 0
    total_images = 0
    avg_loss = 0

    with torch.no_grad():
        for ii, (images, base_labels, type_labels) in enumerate(test_loader):
            type_labels = type_labels.type(torch.LongTensor)
            images = images.type(torch.FloatTensor)
            if gpu_mode:
                images = images.cuda()
                type_labels = type_labels.cuda()

            output_type = transducer_model(images, train_mode=True)
            loss = criterion_type(output_type.contiguous().view(-1, ImageSizeOptions.TOTAL_TYPE_LABELS), type_labels.contiguous().view(-1))

            confusion_matrix.add(output_type.data.contiguous().view(-1, ImageSizeOptions.TOTAL_TYPE_LABELS), type_labels.data.contiguous().view(-1))

            total_loss = loss.item()
            total_images = images.size(0)

            # calculate avg_loss
            avg_loss = total_loss / total_images if total_images else 0

            if (ii + 1) % 100 == 0:
                cm_value = confusion_matrix.value()
                denom = cm_value.sum() if cm_value.sum() > 0 else 1.0

                total_accurate_type = 0
                for i in range(0, ImageSizeOptions.TOTAL_TYPE_LABELS):
                    total_accurate_type = total_accurate_type + cm_value[i][i]

                accuracy_type = (100.0 * total_accurate_type) / denom

                percent_complete = int((100 * (ii+1)) / len(test_loader))
                time_now = time.time()
                mins = int((time_now - start_time) / 60)
                secs = int((time_now - start_time)) % 60

                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO:"
                                 + " BATCH: " + str(ii+1) + "/" + str(len(test_loader))
                                 + " ACCURACY: " + str(round(accuracy_type, 5))
                                 + " LOSS: " + "{:.9f}".format(avg_loss)
                                 + " COMPLETE (" + str(percent_complete) + "%)"
                                 + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                sys.stderr.flush()

    # Confusion matrix
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST LOSS: " + str(avg_loss) + "\n")
    sys.stderr.write("Confusion Matrix:" + "\n")
    sys.stderr.write("            ")
    sys.stderr.flush()
    for label in ImageSizeOptions.decoded_labels:
        sys.stderr.write(str(label) + '    ')
    sys.stderr.write("\n")
    sys.stderr.flush()

    for i, row in enumerate(confusion_matrix.value()):
        sys.stderr.write(str(ImageSizeOptions.decoded_labels[i]) + '   ')
        for j, val in enumerate(row):
            sys.stderr.write("{0:9d}".format(val) + '  ')
        sys.stderr.write("\n")
    sys.stderr.flush()

    cm_value = confusion_matrix.value()
    denom = cm_value.sum() if cm_value.sum() > 0 else 1.0

    total_type_accurate = 0
    for i in range(0, ImageSizeOptions.TOTAL_TYPE_LABELS):
        total_type_accurate = total_type_accurate + cm_value[i][i]

    accuracy = (100.0 * total_type_accurate) / denom

    return {'loss': avg_loss, 'accuracy': accuracy, 'confusion_matrix': confusion_matrix}