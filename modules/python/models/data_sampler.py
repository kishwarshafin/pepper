import torch
from torch.utils.data import sampler
import pandas as pd


class BalancedSampler(sampler.Sampler):
    """
    A sampler that balances the data set.
    Arguments:
        A data set containing images
    """

    def __init__(self, csv_path):
        labels = pd.read_csv(csv_path, header=None)[4]

        self.indices = list(range(len(labels)))
        self.num_samples = len(self.indices)

        label_frequency = labels.value_counts().to_dict()

        data = {'label': labels, 'weight': 0.0}
        weights = pd.DataFrame(data=data)
        weights.loc[weights.label == 0, 'weight'] = 1.0 / label_frequency[0] if 0 in label_frequency.keys() else 0.0
        weights.loc[weights.label == 1, 'weight'] = 1.0 / label_frequency[1] if 1 in label_frequency.keys() else 0.0
        weights.loc[weights.label == 2, 'weight'] = 1.0 / label_frequency[2] if 2 in label_frequency.keys() else 0.0
        weights.loc[weights.label == 3, 'weight'] = 1.0 / label_frequency[3] if 3 in label_frequency.keys() else 0.0
        weights.loc[weights.label == 4, 'weight'] = 1.0 / label_frequency[4] if 4 in label_frequency.keys() else 0.0
        weights.loc[weights.label == 5, 'weight'] = 1.0 / label_frequency[5] if 5 in label_frequency.keys() else 0.0

        weights = weights['weight'].tolist()
        self.weights = torch.DoubleTensor(weights)

    def __iter__(self):
        return (self.indices[i] for i in torch.multinomial(self.weights, self.num_samples))

    def __len__(self):
        return self.num_samples
