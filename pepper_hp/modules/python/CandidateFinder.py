import h5py
import sys
from os.path import isfile, join
from os import listdir
import concurrent.futures
import numpy as np
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
    candidates = sorted(candidates, key=lambda x: (-max(x[10], x[11]), -x[12], -x[6]))

    if len(candidates) > CandidateFinderOptions.MOST_ALLOWED_CANDIDATES_PER_SITE:
        candidates = candidates[0: CandidateFinderOptions.MOST_ALLOWED_CANDIDATES_PER_SITE]

    for i, candidate in enumerate(candidates):
        pos_start, pos_end, ref, alt, alt_type, depth, read_support, \
        read_support_h0, read_support_h1, read_support_h2, alt_prob_h1, alt_prob_h2, non_ref_prob = candidate

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

        if alt_prob_h1 > CandidateFinderOptions.ALT_PROB_THRESHOLD:
            if h1_indx == -1:
                h1_indx = i
                max_h1_prob = alt_prob_h1
            elif max_h1_prob < alt_prob_h1:
                h1_indx = i
                max_h1_prob = alt_prob_h1

        if alt_prob_h2 > CandidateFinderOptions.ALT_PROB_THRESHOLD:
            if h2_indx == -1:
                h2_indx = i
                max_h2_prob = alt_prob_h2
            elif max_h2_prob < alt_prob_h2:
                h2_indx = i
                max_h2_prob = alt_prob_h2

    # for i, candidate in enumerate(candidates):
    #     pos_start, pos_end, ref, alt, alt_type, depth, read_support, \
    #     read_support_h0, read_support_h1, read_support_h2, alt_prob_h1, alt_prob_h2, non_ref_prob = candidate
    #
    #     # make it universal just pick out all
    #     if i in [h1_indx, h2_indx]:
    #         if min_pos_start == -1:
    #             min_pos_start = pos_start
    #         if max_pos_end == -1:
    #             max_pos_end = pos_end
    #         min_pos_start = min(min_pos_start, pos_start)
    #         max_pos_end = max(max_pos_end, pos_end)
    #
    #         if max_pos_end == pos_end:
    #             ref_sequence = ref
    # print(candidates)
    # print(h1_indx, h2_indx)

    selected_alts = []
    selected_dps = []
    selected_alt_prob_h1s = []
    selected_alt_prob_h2s = []
    selected_non_ref_probs = []
    selected_ads = []

    other_alts = []
    other_dps = []
    other_alt_prob_h1s = []
    other_alt_prob_h2s = []
    other_non_ref_probs = []
    other_ads = []
    for i, candidate in enumerate(candidates):
        pos_start, pos_end, ref, alt, alt_type, depth, read_support, \
        read_support_h0, read_support_h1, read_support_h2, alt_prob_h1, alt_prob_h2, non_ref_prob = candidate

        if pos_end < max_pos_end:
            bases_needed = max_pos_end - pos_end
            ref_suffix = ref_sequence[-bases_needed:]
            alt = alt + ref_suffix

        if i in [h1_indx, h2_indx]:
            selected_alts.append(alt)
            selected_dps.append(depth)
            selected_ads.append(read_support)
            selected_alt_prob_h1s.append(alt_prob_h1)
            selected_alt_prob_h2s.append(alt_prob_h2)
            selected_non_ref_probs.append(non_ref_prob)
        else:
            other_alts.append(alt)
            other_dps.append(depth)
            other_ads.append(read_support)
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
    alt_prob_h1s = selected_alt_prob_h1s + other_alt_prob_h1s
    alt_prob_h2s = selected_alt_prob_h2s + other_alt_prob_h2s
    non_ref_probs = selected_non_ref_probs + other_non_ref_probs
    ads = selected_ads + other_ads

    # only report the selected alts
    # alleles = selected_alts
    # dps = selected_dps
    # gts = selected_gts
    # ads = selected_ads
    # print(contig, min_pos_start, max_pos_end, ref_sequence, alleles, genotype, dps, alt_prob_h1s, alt_prob_h2s, ads, overall_non_ref_prob)

    return contig, min_pos_start, max_pos_end, ref_sequence, alleles, genotype, dps, alt_prob_h1s, alt_prob_h2s, non_ref_probs, ads, overall_non_ref_prob


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
        all_predictions_hp1 = []
        all_predictions_hp2 = []

        for i, chunk in enumerate(smaller_chunks):
            with h5py.File(file_name, 'r') as hdf5_file:
                bases_hp1 = hdf5_file['predictions'][contig][chunk_name][chunk]['base_predictions_hp1'][()]
                bases_hp2 = hdf5_file['predictions'][contig][chunk_name][chunk]['base_predictions_hp2'][()]
                positions = hdf5_file['predictions'][contig][chunk_name][chunk]['position'][()]
                indices = hdf5_file['predictions'][contig][chunk_name][chunk]['index'][()]

            if i == 0:
                all_positions = positions
                all_indicies = indices
                all_predictions_hp1 = bases_hp1
                all_predictions_hp2 = bases_hp2
            else:
                all_positions = np.concatenate((all_positions, positions), axis=0)
                all_indicies = np.concatenate((all_indicies, indices), axis=0)
                all_predictions_hp1 = np.concatenate((all_predictions_hp1, bases_hp1), axis=0)
                all_predictions_hp2 = np.concatenate((all_predictions_hp2, bases_hp2), axis=0)

        cpp_candidate_finder = CandidateFinderCPP(contig, contig_start, contig_end)

        # find candidates
        candidate_map = cpp_candidate_finder.find_candidates(bam_file_path,
                                                             reference_file_path,
                                                             contig,
                                                             contig_start,
                                                             contig_end,
                                                             all_positions,
                                                             all_indicies,
                                                             all_predictions_hp1,
                                                             all_predictions_hp2)
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
                                            candidate.alt_prob_h1, candidate.alt_prob_h2, candidate.non_ref_prob))
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
