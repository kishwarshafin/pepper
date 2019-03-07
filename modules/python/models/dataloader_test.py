import os
import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import torchvision.transforms as transforms
import h5py
import torch


class SequenceDataset(Dataset):
    """
    Arguments:
        A CSV file path
    """

    def __init__(self, csv_path, transform=None):
        data_frame = pd.read_csv(csv_path, header=None, dtype=str)
        # assert data_frame[0].apply(lambda x: os.path.isfile(x.split(' ')[0])).all(), \
        #     "Some images referenced in the CSV file were not found"
        self.transform = transforms.Compose([transforms.ToTensor()])
        self.file_info = list(data_frame[0])
        self.file_index = list(data_frame[1])
        self.chromosome_name = list(data_frame[2])

    def __getitem__(self, index):
        # load the image
        chromosome = self.chromosome_name[index]

        hdf5_image = self.file_info[index]
        hdf5_index = int(self.file_index[index])

        hdf5_file = h5py.File(hdf5_image, 'r')

        image_dataset = hdf5_file['image']
        image = image_dataset[hdf5_index]
        image = torch.Tensor(image)

        pos_dataset = hdf5_file['position']
        position = np.array(pos_dataset[hdf5_index], dtype=np.long)

        index_dataset = hdf5_file['index']
        index = np.array(index_dataset[hdf5_index], dtype=np.int16)

        label_dataset = hdf5_file['label']
        label = label_dataset[hdf5_index]
        label = np.array(label, dtype=np.int)

        return image, label, chromosome, position, index

    def __len__(self):
        return len(self.file_info)
