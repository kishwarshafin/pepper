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
        assert data_frame[0].apply(lambda x: os.path.isfile(x.split(' ')[0])).all(), \
            "Some images referenced in the CSV file were not found"
        self.transform = transforms.Compose([transforms.ToTensor()])
        self.file_info = list(data_frame[0])
        self.file_index = list(data_frame[1])

    def __getitem__(self, index):
        # load the image
        hdf5_image = self.file_info[index]
        hdf5_index = int(self.file_index[index])

        hdf5_file = h5py.File(hdf5_image, 'r')
        image_dataset = hdf5_file['image']
        label_dataset = hdf5_file['label']
        image = image_dataset[hdf5_index]
        label = label_dataset[hdf5_index]
        hdf5_file.close()

        image = torch.Tensor(image)
        # label = torch.(image).type(torch.DoubleStorage)
        label = np.array(label, dtype=np.int)
        # label = torch.from_numpy(label).type(torch.LongStorage)

        return image, label

    def __len__(self):
        return len(self.file_info)
