from os.path import isfile, join
from os import listdir
from torch.utils.data import Dataset
from pepper_variant.modules.python.Options import ImageSizeOptions
import concurrent.futures
import torch
import gc
import torchvision.transforms as transforms
import h5py
import gzip
import pickle
import sys
import time
from datetime import datetime
import numpy as np


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-4:] == 'hdf5']
    return file_paths


class SequenceDataset(Dataset):
    """
    Arguments:
        A pkl directory
    """
    def __init__(self, image_directory):
        input_files = get_file_paths_from_directory(image_directory)
        self.all_images = []
        self.all_base_labels = []
        self.all_type_labels = []

        start_time = time.time()
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]" + " INFO: STARTING TO LOAD IMAGES.\n")
        sys.stderr.flush()

        for input_file in input_files:
            with h5py.File(input_file, 'r') as hdf5_file:
                if 'summaries' in hdf5_file:
                    summary_names = list(hdf5_file['summaries'].keys())
                    for summary_name in summary_names:
                        images = hdf5_file['summaries'][summary_name]['images'][()]
                        base_labels = hdf5_file['summaries'][summary_name]['base_labels'][()]
                        type_labels = hdf5_file['summaries'][summary_name]['type_label'][()]
                        self.all_images.extend(images)
                        self.all_base_labels.extend(base_labels)
                        self.all_type_labels.extend(type_labels)

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]" + " INFO: IMAGE LOADING FINISHED.\n")
        sys.stderr.flush()

        time_now = time.time()
        mins = int((time_now - start_time) / 60)
        secs = int((time_now - start_time)) % 60
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]" + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
        sys.stderr.flush()

    def __getitem__(self, index):
        # load the image
        image = np.array(self.all_images[index])
        base_label = self.all_base_labels[index]
        type_label = self.all_type_labels[index]

        return image, base_label, type_label

    def __len__(self):
        return len(self.all_images)


class SequenceDatasetFake(Dataset):
    """
    Arguments:
        A pkl directory
    """
    def __init__(self, image_directory):
        # self.transform = transforms.Compose([transforms.ToTensor()])
        # self.transform = transforms.Compose([])
        input_files = get_file_paths_from_directory(image_directory)
        self.all_contigs = []
        self.all_positions = []
        self.all_depths = []
        self.all_candidates = []
        self.all_candidate_frequency = []
        self.all_images = []
        self.all_base_labels = []
        self.all_type_labels = []

        start_time = time.time()
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]" + " INFO: LOADING: " + str(input_files) + "\n")
        sys.stderr.flush()

        for input_file in input_files:
            with h5py.File(input_file, 'r') as hdf5_file:
                if 'summaries' in hdf5_file:
                    summary_names = list(hdf5_file['summaries'].keys())
                    for summary_name in summary_names:
                        contigs = hdf5_file['summaries'][summary_name]['contigs'][()]
                        positions = hdf5_file['summaries'][summary_name]['positions'][()]
                        depths = hdf5_file['summaries'][summary_name]['depths'][()]
                        candidates = hdf5_file['summaries'][summary_name]['candidates'][()]
                        candidate_frequency = hdf5_file['summaries'][summary_name]['candidate_frequency'][()]
                        images = hdf5_file['summaries'][summary_name]['images'][()]
                        base_labels = hdf5_file['summaries'][summary_name]['base_labels'][()]
                        type_labels = hdf5_file['summaries'][summary_name]['type_label'][()]

                        self.all_contigs.extend(contigs)
                        self.all_positions.extend(positions)
                        self.all_depths.extend(depths)
                        self.all_candidates.extend(candidates)
                        self.all_candidate_frequency.extend(candidate_frequency)
                        self.all_images.extend(images)
                        self.all_base_labels.extend(base_labels)
                        self.all_type_labels.extend(type_labels)

        time_now = time.time()
        mins = int((time_now - start_time) / 60)
        secs = int((time_now - start_time)) % 60
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]" + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
        sys.stderr.flush()

    @staticmethod
    def my_collate(batch):
        contig = [item[0] for item in batch]
        position = [item[1] for item in batch]
        depth = [item[2] for item in batch]
        candidate = [item[3] for item in batch]
        candidate_frequency = [item[4] for item in batch]
        image = [item[5] for item in batch]
        base_label = [item[6] for item in batch]
        type_label = [item[7] for item in batch]

        image = torch.FloatTensor(image)

        return [contig, position, depth, candidate, candidate_frequency, image, base_label, type_label]

    def __getitem__(self, index):
        # load the image
        contig = self.all_contigs[index].decode('UTF-8')

        position = self.all_positions[index]
        depth = self.all_depths[index]
        candidate = self.all_candidates[index]
        candidate_frequency = self.all_candidate_frequency[index]
        image = self.all_images[index]
        base_label = self.all_base_labels[index]
        type_label = self.all_type_labels[index]

        candidate_frequency = [int(x) for x in candidate_frequency]
        candidate = [x for x in candidate]

        base_predictions = np.zeros(ImageSizeOptions.TOTAL_LABELS)
        base_predictions[base_label] = 1

        type_predictions = np.zeros(ImageSizeOptions.TOTAL_TYPE_LABELS)
        type_predictions[type_label] = 1

        # return candidate, image, base_predictions, type_predictions
        return contig, position, depth, candidate, candidate_frequency, image, base_predictions, type_predictions

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