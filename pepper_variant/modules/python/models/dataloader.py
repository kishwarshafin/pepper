from os.path import isfile, join
from os import listdir
from torch.utils.data import Dataset
from pepper_variant.modules.python.Options import ImageSizeOptions
import concurrent.futures
import torch
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
                  and file[-3:] == 'pkl']
    return file_paths


def load_single_file(pickle_files, file_id):
    all_candidates = []
    with open(pickle_files[file_id], "rb") as image_file:
        while True:
            try:
                candidates = pickle.load(image_file)
                all_candidates.extend(candidates)
            except EOFError:
                break

    return all_candidates


def load_all_file(pickle_files):
    all_candidates = []
    total_processes = len(pickle_files)

    with concurrent.futures.ProcessPoolExecutor(max_workers=total_processes) as executor:
        futures = [executor.submit(load_single_file, pickle_files, process_id) for process_id in range(0, total_processes)]

        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                # get the results
                candidates = fut.result()
                all_candidates.extend(candidates)
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    return all_candidates




class SequenceDataset(Dataset):
    """
    Arguments:
        A pkl directory
    """
    def __init__(self, image_directory):
        # self.transform = transforms.Compose([transforms.ToTensor()])
        # self.transform = transforms.Compose([])
        pickle_files = get_file_paths_from_directory(image_directory)
        self.all_candidates = []

        start_time = time.time()
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]" + " INFO: STARTING LOADING PKL.")
        sys.stderr.flush()

        self.all_candidates = load_all_file(pickle_files)

        time_now = time.time()
        mins = int((time_now - start_time) / 60)
        secs = int((time_now - start_time)) % 60
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]" + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
        sys.stderr.flush()
        exit()

    def __getitem__(self, index):
        # load the image
        candidate = self.all_candidates[index]
        image = np.array(candidate.image_matrix)
        type_label = candidate.type_label
        base_label = candidate.base_label

        return image, base_label, type_label

    def __len__(self):
        return len(self.all_candidates)


class SequenceDatasetFake(Dataset):
    """
    Arguments:
        A pkl directory
    """
    def __init__(self, image_directory):
        # self.transform = transforms.Compose([transforms.ToTensor()])
        # self.transform = transforms.Compose([])
        pickle_files = get_file_paths_from_directory(image_directory)
        self.all_candidates = []

        start_time = time.time()
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]" + " INFO: STARTING LOADING PKL.\n")
        sys.stderr.flush()

        self.all_candidates = load_all_file(pickle_files)

        time_now = time.time()
        mins = int((time_now - start_time) / 60)
        secs = int((time_now - start_time)) % 60
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]" + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
        sys.stderr.flush()
        exit()

    @staticmethod
    def my_collate(batch):
        candidate = [item[0] for item in batch]
        images = [item[1] for item in batch]
        base_predictions = [item[2] for item in batch]
        type_predictions = [item[3] for item in batch]
        images = torch.LongTensor(images)

        return [candidate, images, base_predictions, type_predictions]

    def __getitem__(self, index):
        # load the image
        candidate = self.all_candidates[index]
        contig = candidate.contig
        image = np.array(candidate.image_matrix)

        position = candidate.position
        type_label = candidate.type_label
        base_label = candidate.base_label
        candidates = candidate.candidates
        candidate_frequency = candidate.candidate_frequency
        depth = candidate.depth

        base_predictions = np.zeros(ImageSizeOptions.TOTAL_LABELS)
        base_predictions[base_label] = 1

        type_predictions = np.zeros(ImageSizeOptions.TOTAL_TYPE_LABELS)
        type_predictions[type_label] = 1
        # print("######")
        # print(contig)
        #
        # print(position)
        # print(depth)
        # print(candidates)
        # print(candidate_frequency)
        # print(type(image))
        # print(base_predictions)
        # print(type_predictions)
        # print("######")
        # exit(0)

        return candidate, image, base_predictions, type_predictions
        # return contig, depth, candidates, candidate_frequency, image, position, base_predictions, type_predictions

    def __len__(self):
        return len(self.all_candidates)


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