import os
import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import torchvision.transforms as transforms
import h5py
import torch
import pickle


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
        self.dict_info = list(data_frame[2])
        self.record = list(data_frame[3])

    def __getitem__(self, index):
        # load the image
        hdf5_image = self.file_info[index]
        hdf5_index = int(self.file_index[index])
        hdf5_file = h5py.File(hdf5_image, 'r')

        image_dataset = hdf5_file['images']
        image = np.array(image_dataset[hdf5_index], dtype=np.uint8)
        image = self.transform(image)
        image = image.transpose(1, 2)

        dict_path = self.dict_info[index]
        record = self.record[index]

        return image, dict_path, record

    def __len__(self):
        return len(self.file_info)
