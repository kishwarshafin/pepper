import sys
import torch
from tqdm import tqdm
import torchnet.meter as meter
import torch.nn as nn
from torch.utils.data import DataLoader
import numpy as np
from modules.python.models.dataloader import SequenceDataset
from modules.python.TextColor import TextColor
from modules.python.Options import ImageSizeOptions
"""
This script will evaluate a model and return the loss value.

Input:
- A trained model
- A test CSV file to evaluate

Returns:
- Loss value
"""
CLASS_WEIGHTS = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]


def test(data_file, batch_size, gpu_mode, encoder_model, decoder_model, num_workers, gru_layers, hidden_size,
         num_classes=ImageSizeOptions.TOTAL_LABELS):
    # transformations = transforms.Compose([transforms.ToTensor()])

    # data loader
    test_data = SequenceDataset(data_file)
    test_loader = DataLoader(test_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=num_workers,
                             pin_memory=gpu_mode)
    sys.stderr.write(TextColor.CYAN + 'Test data loaded\n')

    # set the evaluation mode of the model
    encoder_model.eval()
    decoder_model.eval()

    class_weights = torch.FloatTensor(CLASS_WEIGHTS)
    # Loss
    test_criterion = nn.CrossEntropyLoss(weight=class_weights)

    if gpu_mode is True:
        test_criterion = test_criterion.cuda()

    # Test the Model
    # sys.stderr.write(TextColor.PURPLE + 'Test starting\n' + TextColor.END)
    confusion_matrix = meter.ConfusionMeter(num_classes)

    total_loss = 0
    total_images = 0
    accuracy = 0

    with torch.no_grad():
        with tqdm(total=len(test_loader), desc='Accuracy: ', leave=True, ncols=100) as pbar:
            for i, (images, labels) in enumerate(test_loader):
                if gpu_mode:
                    # encoder_hidden = encoder_hidden.cuda()
                    images = images.cuda()
                    labels = labels.cuda()

                encoder_hidden = torch.FloatTensor(labels.size(0), gru_layers * 2, hidden_size).zero_()

                if gpu_mode:
                    encoder_hidden = encoder_hidden.cuda()

                # decoder_attentions = torch.zeros(images.size(0), images.size(2), 50)

                loss = 0
                total_seq_length = images.size(2)
                start_index = ImageSizeOptions.CONTEXT_SIZE
                end_index = start_index + (ImageSizeOptions.SEQ_LENGTH - 2 * ImageSizeOptions.CONTEXT_SIZE)

                # from analysis.analyze_png_img import analyze_tensor
                # print(labels[0, :].data.numpy())
                # analyze_tensor(images[0, :, start_index:end_index, :])

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

                    # decoder_attentions[:, seq_index] = attn.view(attn.size(0), -1).data

                    # loss + optimize
                    loss += test_criterion(output_dec, y)
                    confusion_matrix.add(output_dec.data.contiguous().view(-1, num_classes),
                                         y.data.contiguous().view(-1))

                total_loss += loss.item()
                total_images += labels.size(0)

                # print(labels.data.cpu().numpy())
                # plt.matshow(decoder_attentions.numpy()[0])
                # plt.title(str(labels.data.cpu().numpy()[0]), fontsize=12)
                # plt.show()
                # exit()
                # plt.savefig('/data/users/kishwar/train_data/plots/'+"plot_" + str(i) + ".png")
                # exit()

                pbar.update(1)
                cm_value = confusion_matrix.value()
                denom = (cm_value.sum() - cm_value[0][0]) if (cm_value.sum() - cm_value[0][0]) > 0 else 1.0
                accuracy = 100.0 * (cm_value[1][1] + cm_value[2][2] + cm_value[3][3] + cm_value[4][4] +
                                    cm_value[5][5]) / denom
                pbar.set_description("Accuracy: " + str(accuracy))

    avg_loss = total_loss / total_images if total_images else 0

    sys.stderr.write(TextColor.YELLOW+'\nTest Loss: ' + str(avg_loss) + "\n"+TextColor.END)
    sys.stderr.write("Confusion Matrix: \n" + str(confusion_matrix.conf) + "\n" + TextColor.END)

    return {'loss': avg_loss, 'accuracy': accuracy, 'confusion_matrix': str(confusion_matrix.conf)}
