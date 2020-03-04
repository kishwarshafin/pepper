import torch
import torch.nn.functional as F
import torch.nn as nn
from pepper.modules.python.models.resnet import resnet18_custom


def reverse_onehot(one_hot_vector):
    reversed_onehot = one_hot_vector.clone()
    reversed_onehot[one_hot_vector == 0] = 1
    reversed_onehot[one_hot_vector != 0] = 0
    return reversed_onehot


class Attention(nn.Module):
    def __init__(self, dim):
        super(Attention, self).__init__()
        self.linear_out = nn.Linear(dim*2, dim)
        self.mask = None

    def set_mask(self, mask):
        """
        Sets indices to be masked
        Args:
            mask (torch.Tensor): tensor containing indices to be masked
        """
        self.mask = mask

    def forward(self, output, context):
        batch_size = output.size(0)
        hidden_size = output.size(2)
        input_size = context.size(1)

        # (batch, out_len, dim) * (batch, in_len, dim) -> (batch, out_len, in_len)
        attn = torch.bmm(output, context.transpose(1, 2))

        attn = F.softmax(attn.view(-1, input_size), dim=1).view(batch_size, -1, input_size)

        # (batch, out_len, in_len) * (batch, in_len, dim) -> (batch, out_len, dim)
        mix = torch.bmm(attn, context)

        # concat -> (batch, out_len, 2*dim)
        combined = torch.cat((mix, output), dim=2)
        # output -> (batch, out_len, dim)
        output = torch.tanh(self.linear_out(combined.view(-1, 2 * hidden_size))).view(batch_size, -1, hidden_size)

        return output, attn


class EncoderCNN(nn.Module):
    def __init__(self, image_channels):
        """Load the pretrained ResNet-152 and replace top fc layer."""
        super(EncoderCNN, self).__init__()
        self.cnn = resnet18_custom(image_channels)

    def forward(self, images):
        """Extract feature vectors from input images."""
        features = self.cnn(images)

        return features


class EncoderCRNN(nn.Module):
    def __init__(self, image_channels, gru_layers, hidden_size, bidirectional=True):
        super(EncoderCRNN, self).__init__()
        self.cnn_encoder = EncoderCNN(image_channels)
        self.hidden_size = hidden_size
        self.bidirectional = bidirectional
        self.num_layers = gru_layers
        self.gru = nn.GRU(5, hidden_size, num_layers=self.num_layers, bidirectional=bidirectional, batch_first=True)
        self.gru.flatten_parameters()

    def forward(self, x, hidden):
        hidden = hidden.transpose(0, 1).contiguous()

        features_cnn = self.cnn_encoder.forward(x)
        batch_size = features_cnn.size(0)
        seq_len = features_cnn.size(2)
        features_cnn = features_cnn.view(batch_size, seq_len, -1)
        # self.gru.flatten_parameters()
        output_rnn, hidden_rnn = self.gru(features_cnn, hidden)

        if self.bidirectional:
            output_rnn = output_rnn.contiguous()
            output_rnn = output_rnn.view(output_rnn.size(0), output_rnn.size(1), 2, -1) \
                .sum(2).view(output_rnn.size(0), output_rnn.size(1), -1)

        hidden_rnn = hidden_rnn.transpose(0, 1).contiguous()

        return output_rnn, hidden_rnn

    def init_hidden(self, batch_size, num_layers=3, num_directions=2):
        return torch.zeros(batch_size, num_directions * num_layers, self.hidden_size)


class AttnDecoderRNN(nn.Module):
    def __init__(self, hidden_size, gru_layers, num_classes, max_length, seq_len, dropout_p=0.1, bidirectional=True):
        super(AttnDecoderRNN, self).__init__()
        self.hidden_size = hidden_size
        self.num_classes = num_classes
        self.dropout_p = dropout_p
        self.max_length = max_length
        self.bidirectional = bidirectional
        self.embedding = nn.Embedding(self.num_classes, self.hidden_size)
        self.attention = Attention(self.hidden_size)
        self.dropout = nn.Dropout(self.dropout_p)
        self.num_layers = gru_layers
        self.gru = nn.GRU(seq_len, self.hidden_size, num_layers=self.num_layers, batch_first=True,
                          bidirectional=True)
        self.gru.flatten_parameters()
        self.out = nn.Linear(self.hidden_size, self.num_classes)

    def forward_step(self, attention_index_onehot, context_vector, encoder_hidden):
        batch_size = attention_index_onehot.size(0)
        output_gru, hidden_gru = self.gru(attention_index_onehot.view(batch_size, 1, -1), encoder_hidden)
        # print("Attention", attention_index_onehot)

        if self.bidirectional:
            output_gru = output_gru.contiguous()
            output_gru = output_gru.view(output_gru.size(0), output_gru.size(1), 2, -1).sum(2) \
                .view(output_gru.size(0), output_gru.size(1), -1)

        output, attn = self.attention.forward(output_gru, context_vector)

        class_probabilities = self.out(output.contiguous().view(-1, self.hidden_size))

        return class_probabilities, hidden_gru, attn

    def forward(self, attention_index_onehot, context_vector, encoder_hidden):
        encoder_hidden = encoder_hidden.transpose(0, 1).contiguous()
        class_probabilities, hidden, attn = self.forward_step(attention_index_onehot, context_vector,
                                                              encoder_hidden)

        hidden = hidden.transpose(0, 1).contiguous()

        return class_probabilities, hidden, attn
