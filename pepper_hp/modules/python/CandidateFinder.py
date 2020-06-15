import h5py
import time
import sys
from os.path import isfile, join
from os import listdir
import concurrent.futures
import numpy as np
from collections import defaultdict
import operator


BASE_ERROR_RATE = 0.0
label_decoder = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: ''}
label_decoder_ref = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: 'N'}
label_decoder_snp = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: '*'}
MATCH_PENALTY = 4
MISMATCH_PENALTY = 6
GAP_PENALTY = 8
GAP_EXTEND_PENALTY = 2
MIN_SEQUENCE_REQUIRED_FOR_MULTITHREADING = 1
SNP_EVENT = 1
INSERT_EVENT = 2
DELETE_EVENT = 3



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


def group_adjacent_mismatches(mismatches):
    all_groups = []
    current_group = []
    for mismatch in mismatches:
        if len(current_group) == 0:
            current_group.append(mismatch)
        elif abs(current_group[-1][0] - mismatch[0]) == 0:
            # this means the previous position and this position has the same value, meaning it is an insert
            current_group.append(mismatch)
        elif abs(current_group[-1][0] - mismatch[0]) == 1 and mismatch[5] == DELETE_EVENT:
            # this means that this position has a delete that needs to be anchored to the previous position
            current_group.append(mismatch)
        else:
            # otherwise stop the current group and start a new one
            all_groups.append(current_group)
            current_group = [mismatch]

    return all_groups


def mismatch_groups_to_variants(mismatch_group):
    ref_dict = defaultdict()
    allele_dict = defaultdict()
    all_positions = set()

    mismatch_group = sorted(mismatch_group, key=operator.itemgetter(0, 1, 4))
    for pos, indx, ref_base, alt, hp_tag, v_type in mismatch_group:
        all_positions.add((pos, indx))
        ref_dict[(pos, indx)] = ref_base
        allele_dict[(pos, indx, hp_tag)] = alt

    all_positions = list(all_positions)
    all_positions = sorted(all_positions, key=operator.itemgetter(0, 1))
    start_pos = all_positions[0][0]
    start_indx = all_positions[0][1]
    end_pos = all_positions[-1][0]
    start_ref = ref_dict[(start_pos, start_indx)]

    if start_indx > 0 or start_ref == 0:
        # print("GROUP ERROR: ", mismatch_group)
        return None

    ref_allele = []
    alt_allele_1 = []
    alt_allele_2 = []
    for pos, indx in all_positions:
        ref_base = 0
        if indx == 0:
            ref_base = ref_dict[(pos, indx)]
            ref_allele.append(ref_base)

        if (pos, indx, 1) in allele_dict.keys():
            alt_allele_1.append(allele_dict[(pos, indx, 1)])
        else:
            alt_allele_1.append(ref_base)

        if (pos, indx, 2) in allele_dict.keys():
            alt_allele_2.append(allele_dict[(pos, indx, 2)])
        else:
            alt_allele_2.append(ref_base)

    ref_seq = ''.join(label_decoder_ref[i] for i in ref_allele)
    alt1_seq = ''.join(label_decoder[i] for i in alt_allele_1)
    alt2_seq = ''.join(label_decoder[i] for i in alt_allele_2)

    if alt1_seq == ref_seq:
        alt1_seq = ''
    if alt2_seq == ref_seq:
        alt2_seq = ''

    v_type = 'SNP'
    if len(ref_seq) > 1 or max(len(alt1_seq), len(alt2_seq)) > 1:
        v_type = 'INDEL'

    return start_pos, end_pos+1, ref_seq, alt1_seq, alt2_seq, v_type


def get_anchor_positions(base_predictions, ref_seq, indices, positions):
    is_diff = base_predictions != ref_seq
    major_positions = np.where(indices == 0)
    minor_positions = np.where(indices != 0)
    delete_anchors = positions[major_positions][np.where(base_predictions[major_positions] == 0)] - 1
    insert_anchors = positions[minor_positions][np.where(base_predictions[minor_positions] != 0)]

    return list(delete_anchors), list(insert_anchors)

def small_chunk_stitch(input_dir, reference_file_path, contig, small_chunk_keys):
    # for chunk_key in small_chunk_keys:
    all_mismatches = list()
    all_positions = list()
    all_anchors = list()
    position_to_ref_base = defaultdict()

    # Find candidates per chunk
    for file_name, chunk_name in small_chunk_keys:
        all_positions = set()
        highest_index_per_pos = defaultdict(lambda: 0)

        with h5py.File(file_name, 'r') as hdf5_file:
            smaller_chunks = set(hdf5_file['predictions'][contig][chunk_name].keys()) - {'contig_start', 'contig_end'}

        smaller_chunks = sorted(smaller_chunks)

        for chunk in smaller_chunks:
            with h5py.File(file_name, 'r') as hdf5_file:
                bases = hdf5_file['predictions'][contig][chunk_name][chunk]['bases'][()]
                positions = hdf5_file['predictions'][contig][chunk_name][chunk]['position'][()]
                indices = hdf5_file['predictions'][contig][chunk_name][chunk]['index'][()]
                ref_seq = hdf5_file['predictions'][contig][chunk_name][chunk]['ref_seq'][()]

            positions = np.array(positions, dtype=np.int64)
            base_predictions = np.array(bases, dtype=np.int)
            ref_seq = np.array(ref_seq, dtype=np.int)
            delete_anchors, insert_anchors = get_anchor_positions(base_predictions, ref_seq, indices, positions)

            indel_anchors = insert_anchors + delete_anchors
            all_anchors.extend(indel_anchors)

            # as I have carried the reference sequence over, we will get the candidates naturally
            for pos, indx, ref_base, base_pred in zip(positions, indices, ref_seq, base_predictions):
                if indx < 0 or pos < 0:
                    continue

                if (pos, indx) not in all_positions:
                    if ref_base != base_pred:
                        if indx == 0 and ref_base == 0:
                            pass
                        else:
                            all_mismatches.append((pos, indx, ref_base, base_pred))
                            if indx == 0:
                                all_positions.add(pos)

                    if indx == 0:
                        position_to_ref_base[pos] = ref_base

                    all_positions.add((pos, indx))
                    highest_index_per_pos[pos] = max(highest_index_per_pos[pos], indx)

    bad_indices = list()
    # now add all the anchor positions
    for pos in all_anchors:
        if pos not in position_to_ref_base.keys():
            continue
        ref_base = position_to_ref_base[pos]

        if ref_base == 0:
            bad_indices.append(pos)
            continue

        if pos not in all_positions:
            all_mismatches.append((pos, 0, ref_base, ref_base))
            all_positions.add(pos)

    return contig, all_mismatches, bad_indices


def find_candidates(input_dir, reference_file_path, contig, sequence_chunk_keys, threads):
    sequence_chunk_keys = sorted(sequence_chunk_keys)

    all_mismatches = list()
    all_positions = set()

    # generate the dictionary in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        file_chunks = chunks(sequence_chunk_keys, max(2, int(len(sequence_chunk_keys) / threads) + 1))

        futures = [executor.submit(small_chunk_stitch, input_dir, reference_file_path, contig, file_chunk)
                   for file_chunk in file_chunks]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                contig, mismatches, bad_indices = fut.result()
                for pos, indx, ref, alt in mismatches:
                    if (pos, indx) not in all_positions and pos not in bad_indices:
                        all_positions.add((pos, indx))
                        event_type = SNP_EVENT
                        if indx > 0:
                            event_type = INSERT_EVENT
                        elif indx == 0 and alt == 0:
                            event_type = DELETE_EVENT
                        all_mismatches.append((pos, indx, ref, alt, 1, event_type))
            else:
                sys.stderr.write("ERROR IN THREAD: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    all_groups = group_adjacent_mismatches(sorted(all_mismatches, key=operator.itemgetter(0, 1, 4)))

    all_candidate_positions = set()
    candidate_positional_map = defaultdict()

    for mismatch_group in all_groups:
        if not mismatch_group:
            continue

        variant = mismatch_groups_to_variants(mismatch_group)

        if variant is not None:
            _st, _end, ref, alt1, alt2, v_type = variant
        else:
            continue

        if len(alt1) >= 1:
            all_candidate_positions.add(_st)
            candidate_positional_map[_st] = (_st, _end, ref, alt1, v_type)

    return all_candidate_positions, candidate_positional_map
