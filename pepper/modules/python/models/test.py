import sys
import torch
from tqdm import tqdm
import torchnet.meter as meter
import torch.nn as nn
from torch.utils.data import DataLoader
from datetime import datetime
from pepper.modules.python.models.dataloader import SequenceDataset
from pepper.modules.python.Options import ImageSizeOptions, TrainOptions
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
         num_classes=ImageSizeOptions.TOTAL_LABELS, print_details=False):
    # transformations = transforms.Compose([transforms.ToTensor()])

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

    class_weights = torch.Tensor(CLASS_WEIGHTS)
    # Loss
    criterion = nn.CrossEntropyLoss(class_weights)

    if gpu_mode is True:
        criterion = criterion.cuda()

    # Test the Model
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: Test starting\n")
    confusion_matrix = meter.ConfusionMeter(num_classes)

    total_loss = 0
    total_images = 0
    accuracy = 0

    with torch.no_grad():
        with tqdm(total=len(test_loader), desc='Accuracy: ', leave=True, ncols=100) as pbar:
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

                pbar.update(1)
                cm_value = confusion_matrix.value()
                denom = cm_value.sum() if cm_value.sum() > 0 else 1.0
                accuracy = 100.0 * (cm_value[0][0] + cm_value[1][1] + cm_value[2][2]
                                    + cm_value[3][3] + cm_value[4][4]) / denom
                pbar.set_description("Accuracy: " + str(round(accuracy, 5)))

    avg_loss = total_loss / total_images if total_images else 0

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] Test Loss: " + str(avg_loss) + "\n")
    sys.stderr.write("Confusion Matrix: \n" + str(confusion_matrix.conf) + "\n")

    return {'loss': avg_loss, 'accuracy': accuracy, 'confusion_matrix': str(confusion_matrix.conf)}
