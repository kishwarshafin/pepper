from os.path import isfile, join
from os import listdir
from torch.utils.data import Dataset
import torchvision.transforms as transforms
import h5py
import sys
import time
import torch
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

    def __init__(self, image_directory, input_file=None, summary_names=None):
        # self.transform = transforms.Compose([transforms.ToTensor()])
        # self.transform = transforms.Compose([])
        if input_file is None:
            input_files = get_file_paths_from_directory(image_directory)
        else:
            input_files = [input_file]

        self.all_contigs = []
        self.all_positions = []
        self.all_depths = []
        self.all_candidates = []
        self.all_candidate_frequency = []
        self.all_images = []

        for input_file in input_files:
            with h5py.File(input_file, 'r') as hdf5_file:
                if 'summaries' in hdf5_file:
                    if summary_names is None:
                        summary_names = list(hdf5_file['summaries'].keys())
                        for summary_name in summary_names:
                            contigs = hdf5_file['summaries'][summary_name]['contigs'][()]
                            positions = hdf5_file['summaries'][summary_name]['positions'][()]
                            depths = hdf5_file['summaries'][summary_name]['depths'][()]
                            candidates = hdf5_file['summaries'][summary_name]['candidates'][()]
                            candidate_frequency = hdf5_file['summaries'][summary_name]['candidate_frequency'][()]
                            images = hdf5_file['summaries'][summary_name]['images'][()]

                            self.all_contigs.extend(contigs)
                            self.all_positions.extend(positions)
                            self.all_depths.extend(depths)
                            self.all_candidates.extend(candidates)
                            self.all_candidate_frequency.extend(candidate_frequency)
                            self.all_images.extend(images)
                    else:
                        for summary_name in summary_names:
                            contigs = hdf5_file['summaries'][summary_name]['contigs'][()]
                            positions = hdf5_file['summaries'][summary_name]['positions'][()]
                            depths = hdf5_file['summaries'][summary_name]['depths'][()]
                            candidates = hdf5_file['summaries'][summary_name]['candidates'][()]
                            candidate_frequency = hdf5_file['summaries'][summary_name]['candidate_frequency'][()]
                            images = hdf5_file['summaries'][summary_name]['images'][()]

                            self.all_contigs.extend(contigs)
                            self.all_positions.extend(positions)
                            self.all_depths.extend(depths)
                            self.all_candidates.extend(candidates)
                            self.all_candidate_frequency.extend(candidate_frequency)
                            self.all_images.extend(images)


    @staticmethod
    def my_collate(batch):
        contig = [item[0] for item in batch]
        position = [item[1] for item in batch]
        depth = [item[2] for item in batch]
        candidate = [item[3] for item in batch]
        candidate_frequency = [item[4] for item in batch]
        image = [item[5] for item in batch]
        image = np.array(image)
        image = torch.FloatTensor(image)
        # print(type(contig), type(position), type(depth), type(candidate), type(candidate_frequency), type(image))

        return [contig, position, depth, candidate, candidate_frequency, image]

    def __getitem__(self, index):
        # print(index)

        # load the image
        contig = self.all_contigs[index].decode('UTF-8')
        position = self.all_positions[index]
        depth = self.all_depths[index]
        candidate = self.all_candidates[index]
        candidate_frequency = self.all_candidate_frequency[index]
        image = self.all_images[index]


        # return candidate, image, base_predictions, type_predictions
        return contig, position, depth, candidate, candidate_frequency, image

    def __len__(self):
        return len(self.all_images)


class SequenceDatasetHP(Dataset):
    """
    Arguments:
        A HDF5 file path
    """
    def __init__(self, image_directory, file_list=None):
        self.transform = transforms.Compose([transforms.ToTensor()])
        file_image_pair = []

        if file_list is not None:
            hdf_files = file_list
        else:
            hdf_files = get_file_paths_from_directory(image_directory)

        for hdf5_file_path in hdf_files:
            with h5py.File(hdf5_file_path, 'r') as hdf5_file:
                if 'summaries' in hdf5_file:
                    image_names = list(hdf5_file['summaries'].keys())

                    for image_name in image_names:
                        file_image_pair.append((hdf5_file_path, image_name))
                else:
                    sys.stderr.write("WARN: NO IMAGES FOUND IN FILE: " + hdf5_file_path + "\n")

        self.all_images = file_image_pair

    def __getitem__(self, index):
        # load the image
        hdf5_filepath, image_name = self.all_images[index]

        with h5py.File(hdf5_filepath, 'r') as hdf5_file:
            image_hp1 = hdf5_file['summaries'][image_name]['image_hp1'][()]
            image_hp2 = hdf5_file['summaries'][image_name]['image_hp2'][()]

            position = hdf5_file['summaries'][image_name]['position'][()]
            index = hdf5_file['summaries'][image_name]['index'][()]
            contig = hdf5_file['summaries'][image_name]['contig'][()]
            chunk_id = hdf5_file['summaries'][image_name]['chunk_id'][()]
            contig_start = hdf5_file['summaries'][image_name]['region_start'][()]
            contig_end = hdf5_file['summaries'][image_name]['region_end'][()]

        return contig, contig_start, contig_end, chunk_id, image_hp1, image_hp2, position, index

    def __len__(self):
        return len(self.all_images)