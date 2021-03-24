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
        # self.gru_encoder.flatten_parameters()
        # self.gru_decoder.flatten_parameters()
        self.dense1 = nn.Linear(self.hidden_size * 2, self.num_classes)
        # self.dense2 = nn.Linear(self.hidden_size, self.num_classes)

    def forward(self, x, hidden):
        hidden = hidden.transpose(0, 1).contiguous()
        self.gru_encoder.flatten_parameters()
        x_out, hidden_out = self.gru_encoder(x, hidden)
        self.gru_decoder.flatten_parameters()
        x_out, hidden_final = self.gru_decoder(x_out, hidden_out)

        x_out = self.dense1(x_out)
        # x = self.dense2(x)
        # if self.bidirectional:
        #     output_rnn = output_rnn.contiguous()
        #     output_rnn = output_rnn.view(output_rnn.size(0), output_rnn.size(1), 2, -1) \
        #         .sum(2).view(output_rnn.size(0), output_rnn.size(1), -1)

        hidden_final = hidden_final.transpose(0, 1).contiguous()
        return x_out, hidden_final

    def init_hidden(self, batch_size, num_layers, bidirectional=True):
        num_directions = 1
        if bidirectional:
            num_directions = 2

        return torch.zeros(batch_size, num_directions * num_layers, self.hidden_size)


def get_new_gru_model(input_channels, image_features, gru_layers, hidden_size, num_classes=5):
    # get a new model
    transducer_model = TransducerGRU(input_channels, image_features, gru_layers, hidden_size, num_classes,
                                     bidirectional=True)
    return transducer_model


checkpoint = torch.load('/Users/kishwar/Kishwar/software/pepper-private/outputs/PEPPER_SNP_R941_ONT_V4.pkl', map_location='cpu')

hidden_size = checkpoint['hidden_size']
gru_layers = checkpoint['gru_layers']
epochs = checkpoint['epochs']
model_state_dict = checkpoint['model_state_dict']


transducer_model = get_new_gru_model(input_channels=1,
                                     image_features=10,
                                     gru_layers=gru_layers,
                                     hidden_size=hidden_size,
                                     num_classes=15)

from collections import OrderedDict
new_model_state_dict = OrderedDict()

for k, v in model_state_dict.items():
    name = k
    if k[0:7] == 'module.':
        name = k[7:]  # remove `module.`
    new_model_state_dict[name] = v

transducer_model.load_state_dict(new_model_state_dict)
torch.save(transducer_model, '/Users/kishwar/Kishwar/software/pepper-private/outputs/PEPPER_SNP_R941_ONT_V4_TORCH_SAVE.pkl')

x = torch.zeros(1, 100, 10)
h = torch.zeros(1, 2 * 1, 128)

torch.onnx.export(transducer_model, (x, h),
                  '/Users/kishwar/Kishwar/software/pepper-private/outputs/PEPPER_SNP_R941_ONT_V4_TORCH_SAVE' + ".onnx",
                  training=False,
                  opset_version=10,
                  do_constant_folding=True,
                  input_names=['input_image', 'input_hidden'],
                  output_names=['output_pred', 'output_hidden'],
                  dynamic_axes={'input_image': {0: 'batch_size'},
                                'input_hidden': {0: 'batch_size'},
                                'output_pred': {0: 'batch_size'},
                                'output_hidden': {0: 'batch_size'}})