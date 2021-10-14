import torch
import torch.nn as nn
from pepper_variant.modules.python.Options import ImageSizeOptions

class TransducerGRU(nn.Module):
    def __init__(self, image_features, gru_layers, hidden_size, num_classes, num_classes_type, bidirectional=True):
        super(TransducerGRU, self).__init__()
        self.hidden_size = hidden_size
        self.bidirectional = bidirectional
        self.num_layers = gru_layers
        self.num_classes = num_classes
        self.num_classes_type = num_classes_type

        self.linear_1_size = 128
        self.linear_2_size = 64
        self.linear_3_size = 32
        self.linear_4_size = 16
        self.linear_5_size = 8

        # residual block
        self.out_channel = 4
        self.conv1 = nn.Conv2d(1, self.out_channel, kernel_size=3, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(self.out_channel)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = nn.Conv2d(self.out_channel, self.out_channel, kernel_size=3, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(self.out_channel)

        # encoder-only transformer
        self.encoder_layer = nn.TransformerEncoderLayer(d_model=self.out_channel * image_features, nhead=8, batch_first=True)
        self.encoder = nn.TransformerEncoder(self.encoder_layer, num_layers=6)

        self.dropout_1 = nn.Dropout(p=0.1)
        self.dropout_2 = nn.Dropout(p=0.2)
        self.linear_1 = nn.Linear((self.out_channel * image_features) * (ImageSizeOptions.CANDIDATE_WINDOW_SIZE + 1), self.linear_1_size)
        self.linear_2 = nn.Linear(self.linear_1_size, self.linear_2_size)
        self.linear_3 = nn.Linear(self.linear_2_size, self.linear_3_size)
        self.linear_4 = nn.Linear(self.linear_3_size, self.linear_4_size)
        self.linear_5 = nn.Linear(self.linear_4_size, self.linear_5_size)

        self.output_layer = nn.Linear(self.linear_5_size, self.num_classes)

    def forward(self, x, hidden, cell_state, train_mode=False):
        # Reshape for CNN
        x = torch.reshape(x, (x.size()[0], 1, x.size()[1], x.size()[2]))

        # Convolution block
        residual = x
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.conv2(out)
        out = self.bn2(out)
        out += residual
        out = self.relu(out)

        # Reshape for Transformer
        out = torch.reshape(out, (out.size()[0], out.size()[2], out.size()[1] * out.size()[3]))
        # Encoder-only transformer
        out = self.encoder(out)

        # Flatten the output of convolution
        out = torch.flatten(out, start_dim=1, end_dim=2)

        # Linear layer 1
        out = self.linear_1(out)
        out = self.relu(out)
        out = self.dropout_1(out)
        # Linear layer 2
        out = self.linear_2(out)
        out = self.relu(out)
        out = self.dropout_1(out)
        # Linear layer 3
        out = self.linear_3(out)
        out = self.relu(out)
        out = self.dropout_2(out)
        # Linear layer 4
        out = self.linear_4(out)
        out = self.relu(out)
        out = self.dropout_1(out)
        # Linear layer 5
        out = self.linear_5(out)
        out = self.relu(out)
        out = self.output_layer(out)

        if train_mode:
            log_softmax = nn.LogSoftmax(dim=1)
            return log_softmax(out)

        softmax = nn.Softmax(dim=1)
        return softmax(out)

    def init_hidden(self, batch_size, num_layers, bidirectional=True):
        num_directions = 1
        if bidirectional:
            num_directions = 2

        return torch.zeros(batch_size, num_directions * num_layers, self.hidden_size)