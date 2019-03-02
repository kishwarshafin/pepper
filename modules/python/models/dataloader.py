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

    @staticmethod
    def load_dictionary(dictionary_location):
        f = open(dictionary_location, 'rb')
        dict = pickle.load(f)
        f.close()
        return dict

    def __getitem__(self, index):
        # load the image
        hdf5_image = self.file_info[index]
        hdf5_index = int(self.file_index[index])

        hdf5_file = h5py.File(hdf5_image, 'r')
        image = hdf5_file['images'][hdf5_index]
        label = hdf5_file['labels'][hdf5_index]
        hdf5_file.close()

        image = np.array(image, dtype=np.uint8)
        image = self.transform(image)
        image = image.transpose(1, 2)

        label = np.array(label, dtype=np.long)
        label = torch.from_numpy(label).type(torch.LongTensor)

        return image, label

    def __len__(self):
        return len(self.file_info)
