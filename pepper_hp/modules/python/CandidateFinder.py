import h5py
import time
import sys
from os.path import isfile, join
from os import listdir
import concurrent.futures
import numpy as np
from collections import defaultdict
import operator
from pepper_hp.modules.python.Options import CandidateFinderOptions
from pepper_hp.modules.python.CandidateFinderCPP import CandidateFinderCPP


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


def get_index_from_base(base):
    base = base.upper()
    if base == '*':
        return 0
    if base == 'A':
        return 1
    if base == 'C':
        return 2
    if base == 'G':
        return 3
    if base == 'T':
        return 4


def filter_candidate(depth, read_support, read_support_h0, read_support_h1, read_support_h2, alt_prob_h1, alt_prob_h2, non_ref_prob):
    allele_frequency = read_support / max(1.0, depth)
    if alt_prob_h1 >= CandidateFinderOptions.ALT_PROB_THRESHOLD or alt_prob_h2 >= CandidateFinderOptions.ALT_PROB_THRESHOLD:
        return True
    # at this point it means it was not supported by any of the haps, so now look at the non_ref_prob
    if non_ref_prob >= CandidateFinderOptions.NON_REF_PROB_THRESHOLD:
        return True
    # then check the frequency, if too high then definately check
    if allele_frequency >= CandidateFinderOptions.ALLELE_FREQ_THRESHOLD:
        return True
    # now it's extremely unlikely to have this as a true-variant, but still give it one more lifeline
    if allele_frequency >= CandidateFinderOptions.ALLELE_FREQ_THRESHOLD_LAST_RESORT and max(non_ref_prob, alt_prob_h1, alt_prob_h2) >= CandidateFinderOptions.PROB_LAST_RESORT:
        return True

    # otherwise, it's highly unlikely to be a true variant.
    return False


def small_chunk_stitch(reference_file_path, bam_file_path, contig, small_chunk_keys):

    # for chunk_key in small_chunk_keys:
    positional_selected_candidate_list = defaultdict(list)

    # Find candidates per chunk
    for file_name, chunk_name in small_chunk_keys:
        with h5py.File(file_name, 'r') as hdf5_file:
            smaller_chunks = set(hdf5_file['predictions'][contig][chunk_name].keys()) - {'contig_start', 'contig_end'}
            contig_start = hdf5_file['predictions'][contig][chunk_name]['contig_start'][()]
            contig_end = hdf5_file['predictions'][contig][chunk_name]['contig_end'][()]

        cpp_candidate_finder = CandidateFinderCPP(contig, contig_start, contig_end)

        # find candidates
        candidate_map = cpp_candidate_finder.find_candidates(bam_file_path, reference_file_path, contig, contig_start, contig_end)
        smaller_chunks = sorted(smaller_chunks)
        prediction_map_h1 = defaultdict()
        prediction_map_h2 = defaultdict()
        max_index_map = defaultdict()

        for chunk in smaller_chunks:
            with h5py.File(file_name, 'r') as hdf5_file:
                bases = hdf5_file['predictions'][contig][chunk_name][chunk]['base_predictions'][()]
                positions = hdf5_file['predictions'][contig][chunk_name][chunk]['position'][()]
                indices = hdf5_file['predictions'][contig][chunk_name][chunk]['index'][()]
                hp_tag = hdf5_file['predictions'][contig][chunk_name][chunk]['hp_tag'][()]

            positions = np.array(positions, dtype=np.int64)
            base_predictions = np.array(bases, dtype=np.int)

            # as I have carried the reference sequence over, we will get the candidates naturally
            for pos, indx, base_pred in zip(positions, indices, base_predictions):
                if indx < 0 or pos < 0:
                    continue

                if hp_tag == 1:
                    prediction_map_h1[(pos, indx)] = base_pred
                else:
                    prediction_map_h2[(pos, indx)] = base_pred

                if indx == 0:
                    max_index_map[pos] = 1
                else:
                    max_index_map[pos] = max(max_index_map[pos], indx + 1)

        for pos in candidate_map.keys():
            print(pos, chunk_name)
            for candidate in candidate_map[pos]:
                print(candidate.pos_start, candidate.pos_end, candidate.allele.ref, candidate.allele.alt, candidate.allele.alt_type)
                # non-ref prob calcuates the probability of having an alt in that region
                non_ref_prob= 0.0

                alt_prob_h1 = 1.0
                alt_prob_h2 = 1.0
                # see if it's a SNP
                if candidate.allele.alt_type == 1:
                    # see the possibility that this position may
                    pos = candidate.pos_start
                    alt_allele_indx = get_index_from_base(candidate.allele.alt)
                    # calculate probability for HP1
                    prob_alt_h1 = prediction_map_h1[(pos, 0)][alt_allele_indx] / max(1.0, sum(prediction_map_h1[(pos, 0)]))

                    # calculate probability for HP2
                    prob_alt_h2 = prediction_map_h2[(pos, 0)][alt_allele_indx] / max(1.0, sum(prediction_map_h2[(pos, 0)]))

                    # non-ref-prob is the probability that this position can have an alt
                    for indx in range(0, max_index_map[pos]):
                        if indx == 0:
                            ref_allele_indx = get_index_from_base(candidate.allele.ref[0])
                        else:
                            ref_allele_indx = get_index_from_base('*')
                        non_ref_prob_h1 = (sum(prediction_map_h1[(pos, indx)]) - prediction_map_h1[(pos, indx)][ref_allele_indx]) / max(1.0, sum(prediction_map_h1[(pos, indx)]))
                        non_ref_prob_h2 = (sum(prediction_map_h2[(pos, indx)]) - prediction_map_h2[(pos, indx)][ref_allele_indx]) / max(1.0, sum(prediction_map_h2[(pos, indx)]))
                        non_ref_prob = max(non_ref_prob, max(non_ref_prob_h1, non_ref_prob_h2))

                    # calculate probability for HP2
                    alt_prob_h1 = alt_prob_h1 * max(0.01, prob_alt_h1)
                    alt_prob_h2 = alt_prob_h2 * max(0.01, prob_alt_h2)
                # INSERTION
                elif candidate.allele.alt_type == 2:
                    alt_allele = candidate.allele.alt
                    pos = candidate.pos_start
                    indx_lim = min(max_index_map[pos], len(alt_allele))

                    for indx in range(1, indx_lim):
                        alt_allele_indx = get_index_from_base(alt_allele[indx])
                        # print(alt_allele[indx], get_index_from_base(alt_allele[indx]), prediction_map_h1[(pos, indx)], prediction_map_h2[(pos, indx)])
                        prob_alt_h1 = prediction_map_h1[(pos, indx)][alt_allele_indx] / max(1.0, sum(prediction_map_h1[(pos, 0)]))
                        prob_alt_h2 = prediction_map_h2[(pos, indx)][alt_allele_indx] / max(1.0, sum(prediction_map_h2[(pos, 0)]))
                        alt_prob_h1 = alt_prob_h1 * max(0.01, prob_alt_h1)
                        alt_prob_h2 = alt_prob_h2 * max(0.01, prob_alt_h2)

                    # non-ref-prob is the probability that this position can have an alt
                    for indx in range(0, max_index_map[pos]):
                        if indx == 0:
                            ref_allele_indx = get_index_from_base(candidate.allele.ref[0])
                        else:
                            ref_allele_indx = get_index_from_base('*')
                        non_ref_prob_h1 = (sum(prediction_map_h1[(pos, indx)]) - prediction_map_h1[(pos, indx)][ref_allele_indx]) / max(1.0, sum(prediction_map_h1[(pos, indx)]))
                        non_ref_prob_h2 = (sum(prediction_map_h2[(pos, indx)]) - prediction_map_h2[(pos, indx)][ref_allele_indx]) / max(1.0, sum(prediction_map_h2[(pos, indx)]))
                        non_ref_prob = max(non_ref_prob, max(non_ref_prob_h1, non_ref_prob_h2))
                    # alt_probability = prediction_map[(candidate.pos_start)][0]
                # DEL
                elif candidate.allele.alt_type == 3:
                    alt_prob_h1 = 1.0
                    alt_prob_h2 = 1.0
                    for pos in range(candidate.pos_start, candidate.pos_end):
                        ref_allele_indx = get_index_from_base(candidate.allele.ref[pos - candidate.pos_start])
                        non_ref_prob_h1 = (sum(prediction_map_h1[(pos, 0)]) - prediction_map_h1[(pos, 0)][ref_allele_indx]) / max(1.0, sum(prediction_map_h1[(pos, 0)]))
                        non_ref_prob_h2 = (sum(prediction_map_h2[(pos, 0)]) - prediction_map_h2[(pos, 0)][ref_allele_indx]) / max(1.0, sum(prediction_map_h2[(pos, 0)]))
                        non_ref_prob = max(non_ref_prob, max(non_ref_prob_h1, non_ref_prob_h2))

                        if pos > candidate.pos_start:
                            prob_del_h1 = prediction_map_h1[(pos, 0)][0] / max(1.0, sum(prediction_map_h1[(pos, 0)]))
                            prob_del_h2 = prediction_map_h2[(pos, 0)][0] / max(1.0, sum(prediction_map_h2[(pos, 0)]))
                            alt_prob_h1 = alt_prob_h1 * max(0.01, prob_del_h1)
                            alt_prob_h2 = alt_prob_h2 * max(0.01, prob_del_h2)

                if filter_candidate(candidate.depth, candidate.read_support, candidate.read_support_h0, candidate.read_support_h1, candidate.read_support_h2, alt_prob_h1, alt_prob_h2, non_ref_prob):
                    # print(candidate.pos_start, candidate.pos_end, candidate.allele.ref, candidate.allele.alt, candidate.allele.alt_type,
                    #       candidate.depth, candidate.read_support, candidate.read_support_h0, candidate.read_support_h1, candidate.read_support_h2, alt_prob_h1, alt_prob_h2, non_ref_prob)
                    positional_selected_candidate_list[candidate.pos_start].append((candidate.pos_start, candidate.pos_end, candidate.allele.ref, candidate.allele.alt, candidate.allele.alt_type,
                                                                                    candidate.depth, candidate.read_support, candidate.read_support_h0, candidate.read_support_h1, candidate.read_support_h2,
                                                                                    alt_prob_h1, alt_prob_h2, non_ref_prob))

    return positional_selected_candidate_list


def find_candidates(input_dir, reference_file_path, bam_file, contig, sequence_chunk_keys, threads):
    sequence_chunk_keys = sorted(sequence_chunk_keys)
    all_positional_candidates = defaultdict(list)
    all_positions = list()
    # generate the dictionary in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        file_chunks = chunks(sequence_chunk_keys, max(2, int(len(sequence_chunk_keys) / threads) + 1))

        futures = [executor.submit(small_chunk_stitch, reference_file_path, bam_file, contig, file_chunk)
                   for file_chunk in file_chunks]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                positional_candidates = fut.result()

                for pos in positional_candidates.keys():
                    if pos in all_positions:
                        continue
                    all_positional_candidates[pos].extend(positional_candidates[pos])
                    all_positions.append(pos)
            else:
                sys.stderr.write("ERROR IN THREAD: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    return all_positional_candidates
