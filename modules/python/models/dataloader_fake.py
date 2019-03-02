import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import h5py
import torch


class SequenceDataset(Dataset):
    """
    Arguments:
        A CSV file path
    """

    def __init__(self, csv_path, transform=None):
        data_frame = pd.read_csv(csv_path, header=None, dtype=str)
        self.transform = transform

        self.file_info = list(data_frame[0])
        self.file_index = list(data_frame[1])
        self.chromosome_name = list(data_frame[2])

    def __getitem__(self, index):
        chromosome = self.chromosome_name[index]
        hdf5_image = self.file_info[index]
        hdf5_index = int(self.file_index[index])
        hdf5_file = h5py.File(hdf5_image, 'r')

        label_dataset = hdf5_file['label']
        label = np.array(label_dataset[hdf5_index], dtype=np.long)
        label = torch.from_numpy(label).type(torch.LongTensor)

        pos_dataset = hdf5_file['position']
        position = np.array(pos_dataset[hdf5_index], dtype=np.long)

        index_dataset = hdf5_file['index']
        index = np.array(index_dataset[hdf5_index], dtype=np.int16)

        return label, chromosome, position, index

    def __len__(self):
        return len(self.file_info)
