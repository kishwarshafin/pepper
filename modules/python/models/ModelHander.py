import torch
from modules.python.models.Seq2Seq_atn import EncoderCRNN, AttnDecoderRNN


class ModelHandler:
    @staticmethod
    def save_checkpoint(state, filename):
        torch.save(state, filename)

    @staticmethod
    def get_new_model(input_channels, gru_layers, hidden_size, seq_len, num_classes=6):
        # get a new model
        encoder_model = EncoderCRNN(image_channels=input_channels, gru_layers=gru_layers, hidden_size=hidden_size)
        decoder_model = AttnDecoderRNN(hidden_size=hidden_size, gru_layers=gru_layers, num_classes=num_classes,
                                       seq_len=seq_len, max_length=1)
        return encoder_model, decoder_model

    @staticmethod
    def load_optimizer(encoder_optimizer, decoder_optimizer, checkpoint_path, gpu_mode):
        if gpu_mode:
            checkpoint = torch.load(checkpoint_path)
            encoder_optimizer.load_state_dict(checkpoint['encoder_optimizer'])
            for state in encoder_optimizer.state.values():
                for k, v in state.items():
                    if isinstance(v, torch.Tensor):
                        state[k] = v.cuda()

            decoder_optimizer.load_state_dict(checkpoint['decoder_optimizer'])
            for state in decoder_optimizer.state.values():
                for k, v in state.items():
                    if isinstance(v, torch.Tensor):
                        state[k] = v.cuda()
        else:
            checkpoint = torch.load(checkpoint_path, map_location='cpu')
            encoder_optimizer.load_state_dict(checkpoint['encoder_optimizer'])
            decoder_optimizer.load_state_dict(checkpoint['decoder_optimizer'])

        return encoder_optimizer, decoder_optimizer

    @staticmethod
    def load_model_for_training(model_path, input_channels, seq_len, num_classes):
        checkpoint = torch.load(model_path, map_location='cpu')
        hidden_size = checkpoint['hidden_size']
        gru_layers = checkpoint['gru_layers']
        epochs = checkpoint['epochs']
        encoder_model, decoder_model = ModelHandler.get_new_model(input_channels=input_channels,
                                                                  gru_layers=gru_layers,
                                                                  hidden_size=hidden_size,
                                                                  seq_len=seq_len,
                                                                  num_classes=num_classes)
        encoder_state_dict = checkpoint['encoder_state_dict']
        decoder_state_dict = checkpoint['decoder_state_dict']

        from collections import OrderedDict
        new_encoder_state_dict = OrderedDict()
        new_decoder_state_dict = OrderedDict()

        for k, v in encoder_state_dict.items():
            name = k
            if k[0:7] == 'module.':
                name = k[7:]  # remove `module.`
            new_encoder_state_dict[name] = v

        for k, v in decoder_state_dict.items():
            name = k
            if k[0:7] == 'module.':
                name = k[7:]  # remove `module.`
            new_decoder_state_dict[name] = v

        encoder_model.load_state_dict(new_encoder_state_dict)
        decoder_model.load_state_dict(new_decoder_state_dict)
        encoder_model.cpu()
        decoder_model.cpu()

        return encoder_model, decoder_model, hidden_size, gru_layers, epochs

