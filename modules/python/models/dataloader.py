import numpy as np
from os.path import isfile, join
from os import listdir
from torch.utils.data import Dataset
import torchvision.transforms as transforms
import h5py
import torch


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-3:] == 'hdf']
    return file_paths


class SequenceDataset(Dataset):
    """
    Arguments:
        A HDF5 file path
    """
    def __init__(self, image_directory):
        hdf_files = get_file_paths_from_directory(image_directory)
        self.transform = transforms.Compose([transforms.ToTensor()])
        file_image_pair = []
        for hdf5_filepath in hdf_files:
            hdf5_file = h5py.File(hdf5_filepath, 'r')
            image_names = hdf5_file['summaries'].keys()

            for image_name in image_names:
                file_image_pair.append((hdf5_filepath, image_name))

            hdf5_file.close()

        self.all_images = file_image_pair

    def __getitem__(self, index):
        # load the image
        hdf5_filepath, image_name = self.all_images[index]
        hdf5_file = h5py.File(hdf5_filepath, 'r')

        image = hdf5_file['summaries'][image_name]['image']
        label = hdf5_file['summaries'][image_name]['label']
        label = torch.LongTensor(label)
        image = torch.Tensor(image)

        hdf5_file.close()

        return image, label

    def __len__(self):
        return len(self.all_images)
