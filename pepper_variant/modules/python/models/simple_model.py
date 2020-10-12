import torch
import torch.nn as nn


class TransducerGRU(nn.Module):
    def __init__(self, image_channels, image_features, gru_layers, hidden_size, num_classes, bidirectional=True):
        super(TransducerGRU, self).__init__()
        self.hidden_size = hidden_size
        self.bidirectional = bidirectional
        self.num_layers = gru_layers
        self.num_classes = num_classes
        self.gru_encoder = nn.GRU(image_features,
                                  hidden_size,
                                  num_layers=self.num_layers,
                                  bidirectional=bidirectional,
                                  batch_first=True)
        self.gru_decoder = nn.GRU(2 * hidden_size,
                                  hidden_size,
                                  num_layers=self.num_layers,
                                  bidirectional=bidirectional,
                                  batch_first=True)

        self.dense1 = nn.Linear(self.hidden_size * 2, self.num_classes)

    def forward(self, x, hidden):
        hidden = hidden.transpose(0, 1).contiguous()

        self.gru_encoder.flatten_parameters()
        x_out, hidden_out = self.gru_encoder(x, hidden)

        self.gru_decoder.flatten_parameters()
        x_out_final, hidden_final = self.gru_decoder(x_out, hidden_out)

        x_pred = self.dense1(x_out_final)

        hidden_final = hidden_final.transpose(0, 1).contiguous()
        return x_pred, hidden_final

    def init_hidden(self, batch_size, num_layers, bidirectional=True):
        num_directions = 1
        if bidirectional:
            num_directions = 2

        return torch.zeros(batch_size, num_directions * num_layers, self.hidden_size)