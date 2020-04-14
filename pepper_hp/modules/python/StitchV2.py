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


BASE_ERROR_RATE = 0.0
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


def get_confident_positions(alignment):
    cigar_string = alignment.cigar_string.replace('=', 'M').replace('X', 'M')

    cigar_tuples = re.findall(r'(\d+)(\w)', cigar_string)

    grouped_tuples = list()
    prev_len = 0
    prev_op = None
    # group the matches together
    for cigar_len, cigar_op in cigar_tuples:
        if prev_op is None:
            prev_op = cigar_op
            prev_len = int(cigar_len)
        elif prev_op == cigar_op:
            # simply extend the operation
            prev_len += int(cigar_len)
        else:
            grouped_tuples.append((prev_op, prev_len))
            prev_op = cigar_op
            prev_len = int(cigar_len)
    if prev_op is not None:
        grouped_tuples.append((prev_op, prev_len))

    ref_index = alignment.reference_begin
    read_index = 0

    for cigar_op, cigar_len in grouped_tuples:
        if cigar_op == 'M' and cigar_len >= 5:
            return ref_index, read_index

        if cigar_op == 'S':
            read_index += cigar_len
        elif cigar_op == 'I':
            read_index += cigar_len
        elif cigar_op == 'D':
            ref_index += cigar_len
        elif cigar_op == 'M':
            ref_index += cigar_len
            read_index += cigar_len
        else:
            raise ValueError(TextColor.RED + "ERROR: INVALID CIGAR OPERATION ENCOUNTERED WHILTE STITCHING: "
                             + str(cigar_op) + "\n")

    return -1, -1


def alignment_stitch(sequence_chunks):
    sequence_chunks = sorted(sequence_chunks, key=lambda element: (element[1], element[2]))
    contig, running_start, running_end, running_sequence = sequence_chunks[0]
    # if len(running_sequence) < 500:
    #     sys.stderr.write("ERROR: CURRENT SEQUENCE LENGTH TOO SHORT: " + sequence_chunk_keys[0] + "\n")
    #     exit()

    aligner = PEPPER.Aligner(MATCH_PENALTY, MISMATCH_PENALTY, GAP_PENALTY, GAP_EXTEND_PENALTY)
    filter = PEPPER.Filter()
    for i in range(1, len(sequence_chunks)):
        _, this_start, this_end, this_sequence = sequence_chunks[i]
        if this_start < running_end:
            # overlap
            overlap_bases = running_end - this_start
            overlap_bases = overlap_bases + int(overlap_bases * BASE_ERROR_RATE)

            reference_sequence = running_sequence[-overlap_bases:]
            read_sequence = this_sequence[:overlap_bases]

            alignment = PEPPER.Alignment()
            aligner.SetReferenceSequence(reference_sequence, len(reference_sequence))
            aligner.Align_cpp(read_sequence, filter, alignment, 0)

            if alignment.best_score == 0:
                # we are going to say that the left sequence is right
                left_sequence = running_sequence

                # we are going to  say right sequence is also right
                right_sequence = this_sequence

                # but there are 10 'N's as overlaps
                overlap_sequence = 10 * 'N'
                # now append all three parts and we have a contiguous sequence
                running_sequence = left_sequence + overlap_sequence + right_sequence
                running_end = this_end
            else:
                pos_a, pos_b = get_confident_positions(alignment)

                if pos_a == -1 or pos_b == -1:
                    # we are going to say that the left sequence is right
                    left_sequence = running_sequence

                    # we are going to  say right sequence is also right
                    right_sequence = this_sequence

                    # but there are 10 'N's as overlaps
                    overlap_sequence = 10 * 'N'

                    # now append all three parts and we have a contiguous sequence
                    running_sequence = left_sequence + overlap_sequence + right_sequence
                    running_end = this_end
                else:
                    # this is a perfect match so we can simply stitch them
                    # take all of the sequence from the left
                    left_sequence = running_sequence[:-overlap_bases]
                    # get the bases that overlapped
                    overlap_sequence = reference_sequence[:pos_a]
                    # get sequences from current sequence
                    right_sequence = this_sequence[pos_b:]

                    # now append all three parts and we have a contiguous sequence
                    running_sequence = left_sequence + overlap_sequence + right_sequence
                    running_end = this_end
        else:
            # this means there was a gap before this chunk, which could be low read coverage in a small contig.
            running_sequence = running_sequence + this_sequence
            running_end = this_end

    return contig, running_start, running_end, running_sequence


def small_chunk_stitch(file_name, contig, small_chunk_keys):
    # for chunk_key in small_chunk_keys:
    name_sequence_tuples_h1 = list()
    name_sequence_tuples_h2 = list()

    for contig_name, _st, _end in small_chunk_keys:
        chunk_name = contig_name + '-' + str(_st) + '-' + str(_end)

        with h5py.File(file_name, 'r') as hdf5_file:
            contig_start = hdf5_file['predictions'][contig][chunk_name]['contig_start'][()]
            contig_end = hdf5_file['predictions'][contig][chunk_name]['contig_end'][()]

        with h5py.File(file_name, 'r') as hdf5_file:
            smaller_chunks = set(hdf5_file['predictions'][contig][chunk_name].keys()) - {'contig_start', 'contig_end'}

        smaller_chunks = sorted(smaller_chunks)
        all_positions = set()
        base_prediction_dict_h1 = defaultdict()
        base_prediction_dict_h2 = defaultdict()

        for chunk in smaller_chunks:
            with h5py.File(file_name, 'r') as hdf5_file:
                bases_h1 = hdf5_file['predictions'][contig][chunk_name][chunk]['bases_h1'][()]
                bases_h2 = hdf5_file['predictions'][contig][chunk_name][chunk]['bases_h2'][()]
                positions = hdf5_file['predictions'][contig][chunk_name][chunk]['position'][()]
                indices = hdf5_file['predictions'][contig][chunk_name][chunk]['index'][()]

            positions = np.array(positions, dtype=np.int64)
            base_predictions_h1 = np.array(bases_h1, dtype=np.int)
            base_predictions_h2 = np.array(bases_h2, dtype=np.int)

            for pos, indx, base_pred_h1, base_pred_h2 in zip(positions, indices, base_predictions_h1, base_predictions_h2):
                if indx < 0 or pos < 0:
                    continue
                if (pos, indx) not in base_prediction_dict_h1:
                    base_prediction_dict_h1[(pos, indx)] = base_pred_h1
                    base_prediction_dict_h2[(pos, indx)] = base_pred_h2
                    all_positions.add((pos, indx))

        pos_list = sorted(list(all_positions), key=lambda element: (element[0], element[1]))
        dict_fetch = operator.itemgetter(*pos_list)
        predicted_base_labels_h1 = list(dict_fetch(base_prediction_dict_h1))
        predicted_base_labels_h2 = list(dict_fetch(base_prediction_dict_h2))
        sequence_h1 = ''.join([label_decoder[base] for base in predicted_base_labels_h1])
        sequence_h2 = ''.join([label_decoder[base] for base in predicted_base_labels_h2])
        name_sequence_tuples_h1.append((contig, contig_start, contig_end, sequence_h1))
        name_sequence_tuples_h2.append((contig, contig_start, contig_end, sequence_h2))

    name_sequence_tuples_h1 = sorted(name_sequence_tuples_h1, key=lambda element: (element[1], element[2]))
    name_sequence_tuples_h2 = sorted(name_sequence_tuples_h2, key=lambda element: (element[1], element[2]))

    contig, running_start, running_end, running_sequence_h1 = alignment_stitch(name_sequence_tuples_h1)
    contig, running_start, running_end, running_sequence_h2 = alignment_stitch(name_sequence_tuples_h2)
    return contig, running_start, running_end, running_sequence_h1, running_sequence_h2


def create_consensus_sequence(hdf5_file_path, contig, sequence_chunk_keys, threads):
    sequence_chunk_keys = sorted(sequence_chunk_keys)
    sequence_chunk_key_list = list()
    for sequence_chunk_key in sequence_chunk_keys:
        contig, st, end = sequence_chunk_key.split('-')
        sequence_chunk_key_list.append((contig, int(st), int(end)))

    sequence_chunk_key_list = sorted(sequence_chunk_key_list, key=lambda element: (element[1], element[2]))

    sequence_chunks_h1 = list()
    sequence_chunks_h2 = list()
    # generate the dictionary in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        file_chunks = chunks(sequence_chunk_key_list, max(MIN_SEQUENCE_REQUIRED_FOR_MULTITHREADING,
                                                          int(len(sequence_chunk_key_list) / threads) + 1))

        futures = [executor.submit(small_chunk_stitch, hdf5_file_path, contig, file_chunk) for file_chunk in file_chunks]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                contig, contig_start, contig_end, sequence_h1, sequence_h2 = fut.result()
                sequence_chunks_h1.append((contig, contig_start, contig_end, sequence_h1))
                sequence_chunks_h2.append((contig, contig_start, contig_end, sequence_h2))
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    sequence_chunks_h1 = sorted(sequence_chunks_h1, key=lambda element: (element[1], element[2]))
    sequence_chunks_h2 = sorted(sequence_chunks_h2, key=lambda element: (element[1], element[2]))
    contig, contig_start, contig_end, sequence_h1 = alignment_stitch(sequence_chunks_h1)
    contig, contig_start, contig_end, sequence_h2 = alignment_stitch(sequence_chunks_h2)

    return sequence_h1, sequence_h2
