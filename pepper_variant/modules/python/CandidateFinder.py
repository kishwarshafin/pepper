import h5py
import time
import sys
from os.path import isfile, join
from os import listdir
import concurrent.futures
import numpy as np
from collections import defaultdict
import operator
from pepper_variant.modules.python.Options import PEPPERVariantCandidateFinderOptions, ImageSizeOptions
from pepper_variant.modules.python.CandidateFinderCPP import CandidateFinderCPP


BASE_ERROR_RATE = 0.0
label_decoder = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: ''}
label_decoder_ref = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: 'N'}
label_decoder_snp = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: '*'}
MIN_SEQUENCE_REQUIRED_FOR_MULTITHREADING = 1
SNP_EVENT = 1
INSERT_EVENT = 2
DELETE_EVENT = 3


def candidates_to_variants(candidates, contig):
    max_h1_prob = 0.0
    max_h2_prob = 0.0
    h1_indx = -1
    h2_indx = -1
    min_pos_start = -1
    max_pos_end = -1
    ref_sequence = ""
    overall_non_ref_prob = -1.0

    # sort candidates by allele-weight, non-ref prob then allele frequency
    candidates = sorted(candidates, key=lambda x: (-max(x[10], x[11]), -x[6], -x[12]))

    if len(candidates) > PEPPERVariantCandidateFinderOptions.MOST_ALLOWED_CANDIDATES_PER_SITE:
        candidates = candidates[0: PEPPERVariantCandidateFinderOptions.MOST_ALLOWED_CANDIDATES_PER_SITE]

    # found_candidates = False
    for i, candidate in enumerate(candidates):
        pos_start, pos_end, ref, alt, alt_type, depth, read_support, \
        read_support_h0, read_support_h1, read_support_h2, alt_prob, alt_prob_h1, alt_prob_h2, non_ref_prob = candidate

        # valid_candidate = True
        # for base in list(ref):
        #     if base.upper() not in ['A', 'C', 'G', 'T']:
        #         valid_candidate = False
        #         break
        #
        # for base in list(alt):
        #     if base.upper() not in ['A', 'C', 'G', 'T']:
        #         valid_candidate = False
        #         break
        #
        # if not valid_candidate:
        #     continue
        #
        # found_candidates = True

        if overall_non_ref_prob < 0:
            overall_non_ref_prob = non_ref_prob

        overall_non_ref_prob = min(non_ref_prob, overall_non_ref_prob)

        # moving this to a separate loop as we want to do it only on the selected candidates
        if min_pos_start == -1:
            min_pos_start = pos_start
        if max_pos_end == -1:
            max_pos_end = pos_end

        min_pos_start = min(min_pos_start, pos_start)
        max_pos_end = max(max_pos_end, pos_end)

        if max_pos_end == pos_end:
            ref_sequence = ref
        # upto this point

        if alt_prob_h1 > PEPPERVariantCandidateFinderOptions.ALT_PROB_THRESHOLD:
            if h1_indx == -1:
                h1_indx = i
                max_h1_prob = alt_prob_h1
            elif max_h1_prob < alt_prob_h1:
                h1_indx = i
                max_h1_prob = alt_prob_h1

        if alt_prob_h2 > PEPPERVariantCandidateFinderOptions.ALT_PROB_THRESHOLD:
            if h2_indx == -1:
                h2_indx = i
                max_h2_prob = alt_prob_h2
            elif max_h2_prob < alt_prob_h2:
                h2_indx = i
                max_h2_prob = alt_prob_h2

    # if not found_candidates:
    #     return None

    selected_alts = []
    selected_dps = []
    selected_alt_probs = []
    selected_alt_prob_h1s = []
    selected_alt_prob_h2s = []
    selected_non_ref_probs = []
    selected_ads = []

    other_alts = []
    other_dps = []
    other_alt_probs = []
    other_alt_prob_h1s = []
    other_alt_prob_h2s = []
    other_non_ref_probs = []
    other_ads = []
    for i, candidate in enumerate(candidates):
        pos_start, pos_end, ref, alt, alt_type, depth, read_support, \
        read_support_h0, read_support_h1, read_support_h2, alt_prob, alt_prob_h1, alt_prob_h2, non_ref_prob = candidate

        if pos_end < max_pos_end:
            bases_needed = max_pos_end - pos_end
            ref_suffix = ref_sequence[-bases_needed:]
            alt = alt + ref_suffix

        if i in [h1_indx, h2_indx]:
            selected_alts.append(alt)
            selected_dps.append(depth)
            selected_ads.append(read_support)
            selected_alt_probs.append(alt_prob)
            selected_alt_prob_h1s.append(alt_prob_h1)
            selected_alt_prob_h2s.append(alt_prob_h2)
            selected_non_ref_probs.append(non_ref_prob)
        else:
            other_alts.append(alt)
            other_dps.append(depth)
            other_ads.append(read_support)
            other_alt_probs.append(alt_prob)
            other_alt_prob_h1s.append(alt_prob_h1)
            other_alt_prob_h2s.append(alt_prob_h2)
            other_non_ref_probs.append(non_ref_prob)

    indx_list = list()
    for i in [h1_indx, h2_indx]:
        if i > -1:
            indx_list.append(i)

    genotype = [0, 0]
    if len(indx_list) == 1:
        genotype = [0, 1]
    elif len(indx_list) == 2:
        if indx_list[0] == indx_list[1]:
            genotype = [1, 1]
        else:
            genotype = [1, 2]

    alleles = selected_alts + other_alts
    dps = selected_dps + other_dps
    alt_probs = selected_alt_probs + other_alt_probs
    alt_prob_h1s = selected_alt_prob_h1s + other_alt_prob_h1s
    alt_prob_h2s = selected_alt_prob_h2s + other_alt_prob_h2s
    non_ref_probs = selected_non_ref_probs + other_non_ref_probs
    ads = selected_ads + other_ads

    # only report the selected alts
    # alleles = selected_alts
    # dps = selected_dps
    # ads = selected_ads
    # alt_probs = selected_alt_probs
    # alt_prob_h1s = selected_alt_prob_h1s
    # alt_prob_h2s = selected_alt_prob_h2s
    # non_ref_probs = selected_non_ref_probs

    return contig, min_pos_start, max_pos_end, ref_sequence, alleles, genotype, dps, alt_prob_h1s, alt_probs, alt_prob_h2s, non_ref_probs, ads, overall_non_ref_prob



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


def check_alleles(allele):
    allele = allele.upper()
    for base in list(allele):
        if base not in ['A', 'C', 'G', 'T']:
            return False
    return True


def small_chunk_stitch(reference_file_path, bam_file_path, contig, small_chunk_keys):

    # for chunk_key in small_chunk_keys:
    selected_candidate_list = []

    # Find candidates per chunk
    for file_name, chunk_name in small_chunk_keys:
        with h5py.File(file_name, 'r') as hdf5_file:
            smaller_chunks = set(hdf5_file['predictions'][contig][chunk_name].keys()) - {'contig_start', 'contig_end'}
            contig_start = hdf5_file['predictions'][contig][chunk_name]['contig_start'][()]
            contig_end = hdf5_file['predictions'][contig][chunk_name]['contig_end'][()]

        smaller_chunks = sorted(smaller_chunks)

        all_positions = []
        all_indicies = []
        all_predictions = []
        for i, chunk in enumerate(smaller_chunks):
            with h5py.File(file_name, 'r') as hdf5_file:
                bases = hdf5_file['predictions'][contig][chunk_name][chunk]['base_predictions'][()]
                positions = hdf5_file['predictions'][contig][chunk_name][chunk]['position'][()]
                indices = hdf5_file['predictions'][contig][chunk_name][chunk]['index'][()]

            if i == 0:
                all_positions = positions
                all_indicies = indices
                all_predictions = bases
            else:
                all_positions = np.concatenate((all_positions, positions), axis=0)
                all_indicies = np.concatenate((all_indicies, indices), axis=0)
                all_predictions = np.concatenate((all_predictions, bases), axis=0)

        cpp_candidate_finder = CandidateFinderCPP(contig, contig_start, contig_end)

        # find candidates
        candidate_map = cpp_candidate_finder.find_candidates(bam_file_path,
                                                             reference_file_path,
                                                             contig,
                                                             contig_start,
                                                             contig_end,
                                                             all_positions,
                                                             all_indicies,
                                                             all_predictions)

        for pos in candidate_map.keys():
            selected_candidates = []
            found_candidate = False
            for candidate in candidate_map[pos]:
                # print(candidate.pos_start, candidate.pos_end, candidate.allele.ref, candidate.allele.alt, candidate.allele.alt_type,
                #       "(DP", candidate.depth, ", SP: ", candidate.read_support, ", FQ: ", candidate.read_support/candidate.depth, ")",
                #       "H0 supp", candidate.read_support_h0, "H1 supp",  candidate.read_support_h1, "H2 supp",  candidate.read_support_h2,
                #       "ALT prob h1:", candidate.alt_prob_h1, "ALT prob h2:", candidate.alt_prob_h2, "Non-ref prob:", candidate.non_ref_prob)
                found_candidate = True
                selected_candidates.append((candidate.pos_start, candidate.pos_end, candidate.allele.ref, candidate.allele.alt, candidate.allele.alt_type,
                                            candidate.depth, candidate.read_support, candidate.read_support_h0, candidate.read_support_h1, candidate.read_support_h2,
                                            candidate.alt_prob, candidate.alt_prob_h1, candidate.alt_prob_h2, candidate.non_ref_prob))
            if found_candidate:
                variant = candidates_to_variants(list(selected_candidates), contig)
                if variant is not None:
                    selected_candidate_list.append(variant)

    return selected_candidate_list


def find_candidates(input_dir, reference_file_path, bam_file, contig, sequence_chunk_keys, threads):
    sequence_chunk_keys = sorted(sequence_chunk_keys)
    all_selected_candidates = list()
    # generate the dictionary in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        file_chunks = chunks(sequence_chunk_keys, max(2, int(len(sequence_chunk_keys) / threads) + 1))

        futures = [executor.submit(small_chunk_stitch, reference_file_path, bam_file, contig, file_chunk)
                   for file_chunk in file_chunks]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                positional_candidates = fut.result()
                all_selected_candidates.extend(positional_candidates)
            else:
                sys.stderr.write("ERROR IN THREAD: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    # print("TOTAL CANDIDATES IN ", contig, len(all_selected_candidates))
    all_selected_candidates = sorted(all_selected_candidates, key=lambda x: x[1])
    # print("SORTED")
    return all_selected_candidates
