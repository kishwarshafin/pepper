import torch
import torch.nn as nn


class TransducerGRU(nn.Module):
    def __init__(self, image_channels, image_features, gru_layers, hidden_size, num_classes, bidirectional=True):
        super(TransducerGRU, self).__init__()
        self.hidden_size = hidden_size
        self.bidirectional = bidirectional
        self.num_layers = gru_layers
        self.num_classes = num_classes

        self.lstm_1_hidden_size = 128
        self.lstm_2_hidden_size = 256
        self.linear_1_size = 256
        self.linear_2_size = 256
        self.linear_3_size = 256
        self.linear_4_size = 128
        self.linear_5_size = 128

        self.encoder = nn.LSTM(image_features,
                               hidden_size,
                               num_layers=self.num_layers,
                               bidirectional=bidirectional,
                               batch_first=True)
        self.decoder = nn.LSTM(2 * hidden_size,
                               hidden_size,
                               num_layers=self.num_layers,
                               bidirectional=bidirectional,
                               batch_first=True)
        self.dropout_1 = nn.Dropout(p=0.01)
        self.linear_1 = nn.Linear(self.hidden_size * 2, self.linear_1_size)

        self.dropout_2 = nn.Dropout(p=0.1)
        self.linear_2 = nn.Linear(self.linear_1_size, self.linear_2_size)

        self.dropout_3 = nn.Dropout(p=0.1)
        self.linear_3 = nn.Linear(self.linear_2_size, self.linear_3_size)

        self.dropout_4 = nn.Dropout(p=0.2)
        self.linear_4 = nn.Linear(self.linear_3_size, self.linear_4_size)

        self.dropout_5 = nn.Dropout(p=0.1)
        self.linear_5 = nn.Linear(self.linear_4_size, self.linear_5_size)

        self.output_layer = nn.Linear(self.linear_5_size, self.num_classes)

    def forward(self, x, hidden, cell_state, train_mode):
        hidden = hidden.transpose(0, 1).contiguous()
        cell_state = cell_state.transpose(0, 1).contiguous()

        self.encoder.flatten_parameters()
        x, (hidden, cell_state) = self.encoder(x, (hidden, cell_state))

        self.decoder.flatten_parameters()
        x, (hidden, cell_state) = self.decoder(x, (hidden, cell_state))

        x = self.dropout_1(x)
        x = self.linear_1(x)

        x = self.dropout_2(x)
        x = self.linear_2(x)

        x = self.dropout_3(x)
        x = self.linear_3(x)

        x = self.dropout_4(x)
        x = self.linear_4(x)

        x = self.dropout_5(x)
        x = self.linear_5(x)

        x = self.output_layer(x)

        if train_mode:
            log_softmax = nn.LogSoftmax(dim=1)
            return log_softmax(x)

        return x

    def init_hidden(self, batch_size, num_layers, bidirectional=True):
        num_directions = 1
        if bidirectional:
            num_directions = 2

        return torch.zeros(batch_size, num_directions * num_layers, self.hidden_size)