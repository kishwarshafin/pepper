import h5py
import argparse
import sys
from os.path import isfile, join
from os import listdir
import concurrent.futures
import numpy as np
from collections import defaultdict
import operator
from modules.python.TextColor import TextColor
from build import PEPPER
import re


BASE_ERROR_RATE = 1.0
label_decoder = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: ''}
MATCH_PENALTY = 4
MISMATCH_PENALTY = 6
GAP_PENALTY = 8
GAP_EXTEND_PENALTY = 2
MIN_SEQUENCE_REQUIRED_FOR_MULTITHREADING = 2


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path)
                  if isfile(join(directory_path, file)) and file[-2:] == 'h5']
    return file_paths


def chunks(file_names, threads):
    """Yield successive n-sized chunks from l."""
    chunks = []
    for i in range(0, len(file_names), threads):
        chunks.append(file_names[i:i + threads])
    return chunks


def chunks_alignment_sequence(alignment_sequence_pairs, min_length):
    """Yield successive n-sized chunks from l."""
    chunks = []
    for i in range(0, len(alignment_sequence_pairs), min_length):
        chunks.append(alignment_sequence_pairs[i:i + min_length])
    return chunks


def small_chunk_stitch(file_name, contig, small_chunk_keys):
    # for chunk_key in small_chunk_keys:
    all_positions = set()
    base_prediction_dict = defaultdict()

    for contig_name, _st, _end in small_chunk_keys:
        chunk_name = contig_name + '-' + str(_st) + '-' + str(_end)

        with h5py.File(file_name, 'r') as hdf5_file:
            smaller_chunks = set(hdf5_file['predictions'][contig][chunk_name].keys()) - {'contig_start', 'contig_end'}

        smaller_chunks = sorted(smaller_chunks)

        for chunk in smaller_chunks:
            with h5py.File(file_name, 'r') as hdf5_file:
                bases = hdf5_file['predictions'][contig][chunk_name][chunk]['bases'][()]
                positions = hdf5_file['predictions'][contig][chunk_name][chunk]['position'][()]
                indices = hdf5_file['predictions'][contig][chunk_name][chunk]['index'][()]

            positions = np.array(positions, dtype=np.int64)
            base_predictions = np.array(bases, dtype=np.int)

            for pos, indx, base_pred in zip(positions, indices, base_predictions):
                if indx < 0 or pos < 0:
                    continue
                if (pos, indx) not in base_prediction_dict:
                    base_prediction_dict[(pos, indx)] = base_pred
                    all_positions.add((pos, indx))

    return all_positions, base_prediction_dict


def create_consensus_sequence(hdf5_file_path, contig, sequence_chunk_keys, threads):
    sequence_chunk_keys = sorted(sequence_chunk_keys)
    sequence_chunk_key_list = list()
    for sequence_chunk_key in sequence_chunk_keys:
        contig, st, end = sequence_chunk_key.split('-')
        sequence_chunk_key_list.append((contig, int(st), int(end)))

    sequence_chunk_key_list = sorted(sequence_chunk_key_list, key=lambda element: (element[1], element[2]))

    all_positions = set()
    all_base_prediction_dict = defaultdict()
    # generate the dictionary in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        file_chunks = chunks(sequence_chunk_key_list, max(MIN_SEQUENCE_REQUIRED_FOR_MULTITHREADING,
                                                          int(len(sequence_chunk_key_list) / threads) + 1))

        futures = [executor.submit(small_chunk_stitch, hdf5_file_path, contig, file_chunk) for file_chunk in file_chunks]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                positions, base_prediction_dict = fut.result()
                all_positions = all_positions.union(positions)
                all_base_prediction_dict.update(base_prediction_dict)
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    pos_list = sorted(list(all_positions), key=lambda element: (element[0], element[1]))
    dict_fetch = operator.itemgetter(*pos_list)
    predicted_base_labels = list(dict_fetch(all_base_prediction_dict))
    sequence = ''.join([label_decoder[base] for base in predicted_base_labels])

    return sequence
