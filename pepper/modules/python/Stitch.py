import h5py
import sys
from os.path import isfile, join
from os import listdir
import concurrent.futures
import numpy as np
from collections import defaultdict
import operator
from pepper.modules.python.Options import ImageSizeOptions
from datetime import datetime


label_decoder = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: ''}
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


def small_chunk_stitch(contig, small_chunk_keys):
    # for chunk_key in small_chunk_keys:
    base_prediction_dictionary = defaultdict()
    # phred_score_dictionary = defaultdict()
    all_positions = set()
    # ignore first 2 * MIN_IMAGE_OVERLAP bases as they are just overlaps
    buffer_positions = ImageSizeOptions.MIN_IMAGE_OVERLAP * 2

    for file_name, contig_name, _st, _end in small_chunk_keys:
        chunk_name = contig_name + '-' + str(_st) + '-' + str(_end)

        with h5py.File(file_name, 'r') as hdf5_file:
            smaller_chunks = set(hdf5_file['predictions'][contig][chunk_name].keys()) - {'contig_start', 'contig_end'}

        smaller_chunks = sorted(smaller_chunks)

        for i, chunk in enumerate(smaller_chunks):
            with h5py.File(file_name, 'r') as hdf5_file:
                bases = hdf5_file['predictions'][contig][chunk_name][chunk]['bases'][()]
                # phred_score = hdf5_file['predictions'][contig][chunk_name][chunk]['phred_score'][()]
                positions = hdf5_file['predictions'][contig][chunk_name][chunk]['position'][()]
                indices = hdf5_file['predictions'][contig][chunk_name][chunk]['index'][()]

            positions = np.array(positions, dtype=np.int64)
            base_predictions = np.array(bases, dtype=np.int)
            # if _st == 107900:
            #     print(positions)

            for pos, indx, base_pred in zip(positions, indices, base_predictions):
                # not take the first buffer bases for every chunk that has an overlap to the last chunk
                if _st > 0 and pos <= _st + buffer_positions:
                    continue

                if indx < 0 or pos < 0:
                    continue

                base_prediction_dictionary[(pos, indx)] = base_pred
                # phred_score_dictionary[(pos, indx)] = base_score + 33
                all_positions.add((pos, indx))

    if len(all_positions) == 0:
        return -1, -1, ''

    pos_list = sorted(list(all_positions), key=lambda element: (element[0], element[1]))
    dict_fetch = operator.itemgetter(*pos_list)

    # weird but python has no casting between np.int64 to list
    if len(pos_list) > 1:
        predicted_base_labels = list(dict_fetch(base_prediction_dictionary))
        # predicted_phred_scores = list(dict_fetch(phred_score_dictionary))
    else:
        predicted_base_labels = [dict_fetch(base_prediction_dictionary)]
        # predicted_phred_scores = [dict_fetch(phred_score_dictionary)]

    sequence = ''.join([label_decoder[base] for base in predicted_base_labels])
    # phred_score = ''.join([chr(phred_score) for base, phred_score in zip(predicted_base_labels, predicted_phred_scores) if base != 0])
    first_pos = pos_list[0][0]
    last_pos = pos_list[-1][0]
    return first_pos, last_pos, sequence


def create_consensus_sequence(contig, sequence_chunk_keys, threads):
    sequence_chunk_keys = sorted(sequence_chunk_keys, key=lambda element: (element[1]))

    sequence_chunk_key_list = list()
    for file_name, sequence_chunk_key, contig_start, contig_end in sequence_chunk_keys:
        sequence_chunk_key_list.append((file_name, contig, int(contig_start), int(contig_end)))

    sequence_chunk_key_list = sorted(sequence_chunk_key_list, key=lambda element: (element[2], element[3]))

    sequence_chunks = list()

    # generate the dictionary in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        file_chunks = chunks(sequence_chunk_key_list, max(MIN_SEQUENCE_REQUIRED_FOR_MULTITHREADING,
                                                          int(len(sequence_chunk_key_list) / threads) + 1))

        futures = [executor.submit(small_chunk_stitch, contig, chunk) for chunk in file_chunks]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                first_pos, last_pos, sequence = fut.result()
                if first_pos != -1 and last_pos != -1:
                    sequence_chunks.append((first_pos, last_pos, sequence))
            else:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]  ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    sequence_chunks = sorted(sequence_chunks, key=lambda element: (element[0], element[1]))
    stitched_sequence = ''
    for first_pos, last_pos, sequence in sequence_chunks:
        stitched_sequence = stitched_sequence + sequence

    return stitched_sequence
