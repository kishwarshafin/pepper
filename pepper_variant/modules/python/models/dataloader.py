from os.path import isfile, join
from os import listdir
from torch.utils.data import Dataset
from pepper_variant.modules.python.Options import ImageSizeOptions
import torchvision.transforms as transforms
from collections import defaultdict
import pandas as pd
import h5py
import sys
import numpy as np


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
    def __init__(self, csv_file):
        self.transform = transforms.Compose([transforms.ToTensor()])
        self.hdf_filenames = defaultdict()
        self.region_names = defaultdict()

        print(csv_file)
        csv_file = pd.read_csv(csv_file, index_col=0, header=None)

        self.all_images = csv_file
        print(type(self.all_images))
        exit()

    def __getitem__(self, index):
        # load the image
        hdf5_filepath, region_name, indx = self.all_images[index]

        with h5py.File(hdf5_filepath, 'r') as hdf5_file:
            image = hdf5_file['summaries'][region_name]['images'][indx][()]
            type_label = hdf5_file['summaries'][region_name]['type_labels'][indx][()]
            base_label = hdf5_file['summaries'][region_name]['base_labels'][indx][()]

        return image, base_label, type_label

    def __len__(self):
        index = self.all_images.index
        number_of_rows = len(index)
        return len(number_of_rows)


class SequenceDatasetFake(Dataset):
    """
    Arguments:
        A HDF5 file path
    """
    def __init__(self, image_directory):
        self.transform = transforms.Compose([transforms.ToTensor()])
        file_image_pair = []

        hdf_files = get_file_paths_from_directory(image_directory)

        for hdf5_file_path in hdf_files:
            with h5py.File(hdf5_file_path, 'r') as hdf5_file:
                if 'summaries' in hdf5_file:
                    region_names = list(hdf5_file['summaries'].keys())
                    for region_name in region_names:
                        image_shape = hdf5_file['summaries'][region_name]['images'].shape[0]

                        for index in range(0, image_shape):
                            file_image_pair.append((hdf5_file_path, region_name, index))
                else:
                    sys.stderr.write("WARN: NO IMAGES FOUND IN FILE: " + hdf5_file_path + "\n")

        self.all_images = file_image_pair

    def __getitem__(self, index):
        # load the image
        hdf5_filepath, region_name, index = self.all_images[index]

        with h5py.File(hdf5_filepath, 'r') as hdf5_file:
            image = hdf5_file['summaries'][region_name]['images'][index][()]
            position = hdf5_file['summaries'][region_name]['positions'][index][()]
            type_label = hdf5_file['summaries'][region_name]['type_labels'][index][()]
            base_label = hdf5_file['summaries'][region_name]['base_labels'][index][()]

            contig = hdf5_file['summaries'][region_name]['contig'][()]
            region_start = hdf5_file['summaries'][region_name]['region_start'][()]
            region_stop = hdf5_file['summaries'][region_name]['region_end'][()]

            base_predictions = np.zeros((base_label.size, ImageSizeOptions.TOTAL_LABELS))
            base_predictions[np.arange(base_label.size), base_label] = 1

            type_predictions = np.zeros((type_label.size, ImageSizeOptions.TOTAL_TYPE_LABELS))
            type_predictions[np.arange(type_label.size), type_label] = 1

        return contig, region_start, region_stop, image, position, base_predictions, type_predictions

    def __len__(self):
        return len(self.all_images)


class SequenceDatasetHP(Dataset):
    """
    Arguments:
        A HDF5 file path
    """
    def __init__(self, image_directory):
        self.transform = transforms.Compose([transforms.ToTensor()])
        file_image_pair = []

        hdf_files = get_file_paths_from_directory(image_directory)

        for hdf5_file_path in hdf_files:
            with h5py.File(hdf5_file_path, 'r') as hdf5_file:
                if 'summaries' in hdf5_file:
                    image_names = list(hdf5_file['summaries'].keys())

                    for image_name in image_names:
                        # for hp_tag 1
                        file_image_pair.append((hdf5_file_path, image_name, 1))
                        # for hp_tag 2
                        file_image_pair.append((hdf5_file_path, image_name, 2))
                else:
                    sys.stderr.write("WARN: NO IMAGES FOUND IN FILE: " + hdf5_file_path + "\n")

        self.all_images = file_image_pair

    def __getitem__(self, index):
        # load the image
        hdf5_filepath, image_name, hp_tag = self.all_images[index]

        with h5py.File(hdf5_filepath, 'r') as hdf5_file:
            if hp_tag == 1:
                image = hdf5_file['summaries'][image_name]['image_hp1'][()]
                label = hdf5_file['summaries'][image_name]['label_hp1'][()]
            else:
                image = hdf5_file['summaries'][image_name]['image_hp2'][()]
                label = hdf5_file['summaries'][image_name]['label_hp2'][()]

        return image, label

    def __len__(self):
        return len(self.all_images)