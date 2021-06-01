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
    all_images = []
    all_type_labels = []
    all_base_labels = []
    hdf_files = get_file_paths_from_directory(input_directory)
    record_index = 0
    print("HDF5 FILES: ", hdf_files)
    for i, hdf5_file_path in enumerate(hdf_files):
        print("PROCESSING: ", i + 1, "/", len(hdf_files))
        with h5py.File(hdf5_file_path, 'r') as hdf5_file:
            if 'summaries' in hdf5_file:
                region_names = list(hdf5_file['summaries'].keys())

                for region_name in region_names:
                    if record_index == 0:
                        all_images = hdf5_file['summaries'][region_name]['images'][()]
                        all_type_labels = hdf5_file['summaries'][region_name]['type_labels'][()]
                        all_base_label = hdf5_file['summaries'][region_name]['base_labels'][()]
                    else:
                        all_images = np.append(all_images, hdf5_file['summaries'][region_name]['images'][()], axis=0)
                        all_type_labels = np.append(all_type_labels, hdf5_file['summaries'][region_name]['type_labels'][()], axis=0)
                        all_base_label = np.append(all_base_label, hdf5_file['summaries'][region_name]['base_labels'][()], axis=0)
                    record_index += 1

        print("IMAGE SHAPE:", all_images.shape)
        print("BASE SHAPE:", all_base_label.shape)
        print("TYPE SHAPE:", all_type_labels.shape)

    output_hdf5_filename = output_directory + "/" + "Merged_file.hdf"
    with h5py.File(output_hdf5_filename, 'w') as output_hdf5_file:
        output_hdf5_file.create_dataset("images", all_images.shape, dtype=np.uint8, compression="gzip", data=all_images)
        output_hdf5_file.create_dataset("base_labels", all_base_label.shape, dtype=np.uint8, data=all_base_label)
        output_hdf5_file.create_dataset("type_labels", all_type_labels.shape, dtype=np.uint8, data=all_type_labels)

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

