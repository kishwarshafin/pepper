import argparse
from os.path import isfile, join
from os import listdir
from torch.utils.data import Dataset
from pepper_variant.modules.python.Options import ImageSizeOptions
import torchvision.transforms as transforms
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


def merge_hdf5_files(input_directory, output_directory):
    hdf_files = get_file_paths_from_directory(input_directory)
    record_index = 0
    output_hdf5_filename = output_directory + "/" + "Merged_file.hdf"
    output_hdf5_file = h5py.File(output_hdf5_filename, 'w')

    print("HDF5 FILES: ", hdf_files)
    for i, hdf5_file_path in enumerate(hdf_files):
        print("PROCESSING: ", i + 1, "/", len(hdf_files))
        with h5py.File(hdf5_file_path, 'r') as hdf5_file:
            if 'summaries' in hdf5_file:
                region_names = list(hdf5_file['summaries'].keys())

                for region_name in region_names:
                    image_shape = hdf5_file['summaries'][region_name]['images'].shape[0]

                    for image_index in range(0, image_shape):
                        output_hdf5_file[str(record_index)] = hdf5_file_path + "," + region_name + "," + str(image_index)

                        record_index += 1

        print("TOTAL RECORDS:", record_index)

    print("FILES SAVED.")


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--image_directory",
        "-i",
        type=str,
        required=True,
        help="Path to a directory containing hdf5 files."
    )
    parser.add_argument(
        "--output_directory",
        "-o",
        type=str,
        required=True,
        help="Path to a output directory."
    )
    FLAGS, unparsed = parser.parse_known_args()

    merge_hdf5_files(FLAGS.image_directory, FLAGS.output_directory)

