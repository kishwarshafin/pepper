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

        self.lstm_1_hidden_size = 256
        self.lstm_2_hidden_size = 256
        self.linear_1_size = 256
        self.linear_2_size = 256
        self.linear_3_size = 256
        self.linear_4_size = 128
        self.linear_5_size = 128

        self.encoder = nn.LSTM(image_features,
                               self.lstm_1_hidden_size,
                               num_layers=self.num_layers,
                               bidirectional=bidirectional,
                               batch_first=True)
        self.decoder = nn.LSTM(2 * self.lstm_2_hidden_size,
                               self.lstm_2_hidden_size,
                               num_layers=self.num_layers,
                               bidirectional=bidirectional,
                               batch_first=True)
        self.dropout_1 = nn.Dropout(p=0.1)

        self.linear_1 = nn.Linear((self.lstm_2_hidden_size * 2) * (ImageSizeOptions.CANDIDATE_WINDOW_SIZE + 1), self.linear_1_size)
        self.relu_1 = nn.ReLU()

        self.dropout_2 = nn.Dropout(p=0.2)
        self.linear_2 = nn.Linear(self.linear_1_size, self.linear_2_size)
        self.relu_2 = nn.ReLU()

        self.dropout_3 = nn.Dropout(p=0.1)
        self.linear_3 = nn.Linear(self.linear_2_size, self.linear_3_size)
        self.relu_3 = nn.ReLU()

        self.dropout_4 = nn.Dropout(p=0.2)
        self.linear_4 = nn.Linear(self.linear_3_size, self.linear_4_size)
        self.relu_4 = nn.ReLU()

        # self.dropout_4_type = nn.Dropout(p=0.2)
        # self.linear_4_type = nn.Linear(self.linear_3_size, self.linear_4_size)
        # self.relu_4_type = nn.ReLU()

        self.dropout_5 = nn.Dropout(p=0.1)
        self.linear_5 = nn.Linear(self.linear_4_size, self.linear_5_size)
        self.relu_5 = nn.ReLU()

        # self.dropout_5_type = nn.Dropout(p=0.1)
        # self.linear_5_type = nn.Linear(self.linear_4_size, self.linear_5_size)
        # self.relu_5_type = nn.ReLU()

        self.output_layer = nn.Linear(self.linear_5_size, self.num_classes)
        # self.output_layer_type = nn.Linear(self.linear_5_size, self.num_classes_type)

    def forward(self, x, hidden, cell_state, train_mode=False):
        hidden = hidden.transpose(0, 1).contiguous()
        cell_state = cell_state.transpose(0, 1).contiguous()

        self.encoder.flatten_parameters()
        x, (hidden, cell_state) = self.encoder(x, (hidden, cell_state))

        self.decoder.flatten_parameters()
        x, (hidden, cell_state) = self.decoder(x, (hidden, cell_state))
        x = self.dropout_1(x)
        x = torch.flatten(x, start_dim=1, end_dim=2)

        x = self.linear_1(x)
        x = self.relu_1(x)
        x = self.dropout_2(x)

        x = self.linear_2(x)
        x = self.relu_2(x)
        x = self.dropout_3(x)

        x = self.linear_3(x)
        x = self.relu_3(x)
        x = self.dropout_4(x)

        x_base = self.linear_4(x)
        x_base = self.relu_4(x_base)
        x_base = self.dropout_5(x_base)

        x_base = self.linear_5(x_base)
        x_base = self.relu_5(x_base)
        x_base = self.output_layer(x_base)

        if train_mode:
            log_softmax = nn.LogSoftmax(dim=1)
            return log_softmax(x_base)

        softmax = nn.Softmax(dim=1)
        return softmax(x_base)

        # x_type = self.linear_4_type(x)
        # x_type = self.relu_4_type(x_type)
        # x_type = self.dropout_5_type(x_type)
        #
        # x_type = self.linear_5_type(x_type)
        # x_type = self.relu_5_type(x_type)
        # x_type = self.output_layer_type(x_type)

        # if train_mode:
        #     log_softmax = nn.LogSoftmax(dim=1)
        #     return log_softmax(x_base), log_softmax(x_type)
        #
        # softmax = nn.Softmax(dim=1)
        # return softmax(x_base), softmax(x_type)

    def init_hidden(self, batch_size, num_layers, bidirectional=True):
        num_directions = 1
        if bidirectional:
            num_directions = 2

        return torch.zeros(batch_size, num_directions * num_layers, self.hidden_size)