import torch
import torch.nn.functional as F
import torch.nn as nn


class TransducerGRU(nn.Module):
    def __init__(self, image_channels, image_features, gru_layers, hidden_size, num_classes, bidirectional=True):
        super(TransducerGRU, self).__init__()
        self.hidden_size = hidden_size
        self.bidirectional = bidirectional
        self.num_layers = gru_layers
        self.num_classes = num_classes
        self.gru = nn.GRU(image_features,
                          hidden_size,
                          num_layers=self.num_layers,
                          bidirectional=bidirectional,
                          batch_first=True)
        self.gru.flatten_parameters()
        self.dense1 = nn.Linear(self.hidden_size * 2, self.num_classes)
        # self.dense2 = nn.Linear(self.hidden_size, self.num_classes)

    def forward(self, x, hidden):
        hidden = hidden.transpose(0, 1).contiguous()
        # self.gru.flatten_parameters()
        x, hidden = self.gru(x, hidden)
        x = self.dense1(x)
        # x = self.dense2(x)
        # if self.bidirectional:
        #     output_rnn = output_rnn.contiguous()
        #     output_rnn = output_rnn.view(output_rnn.size(0), output_rnn.size(1), 2, -1) \
        #         .sum(2).view(output_rnn.size(0), output_rnn.size(1), -1)

        # hidden_rnn = hidden_rnn.transpose(0, 1).contiguous()
        return x

    def init_hidden(self, batch_size, num_layers=3, num_directions=2):
        return torch.zeros(batch_size, num_directions * num_layers, self.hidden_size)