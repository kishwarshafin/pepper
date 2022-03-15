import h5py
import sys
import math
from os.path import isfile, join
from os import listdir
import concurrent.futures
import numpy as np
from collections import defaultdict
from pepper_variant.modules.python.Options import PEPPERVariantCandidateFinderOptions, ImageSizeOptions
from pepper_variant.build import PEPPER_VARIANT


BASE_ERROR_RATE = 0.0
label_decoder = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: ''}
label_decoder_ref = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: 'N'}
label_decoder_snp = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: '*'}
MIN_SEQUENCE_REQUIRED_FOR_MULTITHREADING = 1
SNP_EVENT = 1
INSERT_EVENT = 2
DELETE_EVENT = 3


def candidates_to_variants(candidates, contig, freq_based, freq):
    max_h1_prob = 0.0
    max_h2_prob = 0.0
    h1_indx = -1
    h2_indx = -1
    min_pos_start = -1
    max_pos_end = -1
    ref_sequence = ""
    overall_non_ref_prob = -1.0

    # sort candidates by allele-weight, non-ref prob then allele frequency
    candidates = sorted(candidates, key=lambda x: (-max(x[7], x[8]), -x[9], -x[6]))

    if len(candidates) > PEPPERVariantCandidateFinderOptions.MOST_ALLOWED_CANDIDATES_PER_SITE:
        if not freq_based:
            candidates = candidates[0: PEPPERVariantCandidateFinderOptions.MOST_ALLOWED_CANDIDATES_PER_SITE]

    # found_candidates = False
    for i, candidate in enumerate(candidates):
        pos_start, pos_end, ref, alt, alt_type, depth, read_support, alt_prob_h1, alt_prob_h2, non_ref_prob = candidate

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
        pos_start, pos_end, ref, alt, alt_type, depth, read_support, alt_prob_h1, alt_prob_h2, non_ref_prob = candidate

        if pos_end < max_pos_end:
            bases_needed = max_pos_end - pos_end
            ref_suffix = ref_sequence[-bases_needed:]
            alt = alt + ref_suffix

        if i in [h1_indx, h2_indx]:
            selected_alts.append(alt)
            selected_dps.append(depth)
            selected_ads.append(read_support)
            selected_alt_probs.append(max(alt_prob_h1, alt_prob_h2))
            selected_alt_prob_h1s.append(alt_prob_h1)
            selected_alt_prob_h2s.append(alt_prob_h2)
            selected_non_ref_probs.append(non_ref_prob)
        else:
            other_alts.append(alt)
            other_dps.append(depth)
            other_ads.append(read_support)
            other_alt_probs.append(max(alt_prob_h1, alt_prob_h2))
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

    if freq_based:
        alleles = selected_alts + other_alts
        dps = selected_dps + other_dps
        alt_probs = selected_alt_probs + other_alt_probs
        alt_prob_h1s = selected_alt_prob_h1s + other_alt_prob_h1s
        alt_prob_h2s = selected_alt_prob_h2s + other_alt_prob_h2s
        non_ref_probs = selected_non_ref_probs + other_non_ref_probs
        ads = selected_ads + other_ads
    else:
        # only report the selected alts
        alleles = selected_alts
        dps = selected_dps
        ads = selected_ads
        alt_probs = selected_alt_probs
        alt_prob_h1s = selected_alt_prob_h1s
        alt_prob_h2s = selected_alt_prob_h2s
        non_ref_probs = selected_non_ref_probs

    return contig, min_pos_start, max_pos_end, ref_sequence, alleles, genotype, dps, alt_probs, alt_prob_h1s, alt_prob_h2s, non_ref_probs, ads, overall_non_ref_prob


def candidates_to_variants_snp(candidates, contig, freq_based, freq):
    ref_sequence = ""

    min_pos_start, max_pos_end, genotype = 0, 0, 0
    reference_sequence = ""
    alleles, dps, ads = [], [], []
    genotype = 0
    allele_probability, genotype_probability = 0, 0

    for i, candidate in enumerate(candidates):
        pos_start, pos_end, ref, alt, alt_type, depth, read_support, gt, allele_prob, genotype_prob = candidate
        if i == 0:
            min_pos_start = pos_start
            max_pos_end = pos_end
            reference_sequence = ref
            alleles.append(alt)
            dps.append(depth)
            ads.append(read_support)
            genotype = gt
            allele_probability = allele_prob
            genotype_probability = genotype_prob
        else:
            alleles.append(alt)
            dps.append(depth)
            ads.append(read_support)

            if min_pos_start != pos_start:
                print("MIN POS START DIDN'T MATCH: ", candidate)
            if max_pos_end != pos_end:
                print("MIN POS START DIDN'T MATCH: ", candidate)
            if ref != reference_sequence:
                print("REFERENCE DIDN'T MATCH: ", candidate)
            if gt != genotype:
                print("GENOTYPE DIDN'T MATCH: ", candidate)
            if allele_probability != allele_prob:
                print("ALLELE PROBABILITY DIDN'T MATCH", candidate)
            if genotype_prob != genotype_prob:
                print("GENOTYPE PROBABILITY DIDN'T MATCH", candidate)

    final_genotype = [0, 0]
    if genotype == 1:
        if len(alleles) == 1:
            final_genotype = [0, 1]
        elif len(alleles) > 1:
            final_genotype = [1, 2]
    elif genotype == 2:
        if len(alleles) == 1:
            final_genotype = [1, 1]
        elif len(alleles) > 1:
            final_genotype = [1, 2]

    return contig, min_pos_start, max_pos_end, reference_sequence, alleles, dps, ads, final_genotype, allele_probability, genotype_probability


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


def get_genotype_from_base(ref_base, base_prediction1, base_prediction2):
    if base_prediction1 == 'R':
        base_prediction1 = ref_base

    if base_prediction2 == 'R':
        base_prediction2 = ref_base

    if ref_base == base_prediction1 or ref_base == base_prediction2:
        if base_prediction1 == base_prediction2:
            return [0, 0]
        else:
            return [0, 1]
    elif base_prediction1 == base_prediction2:
        return [1, 1]
    else:
        return [0, 1]


def repeat_annotation(sequence, kmer_size):
    """
    Annotate sequence with k-mer repeats.
    """
    max_observed_repeats = [1 for i in range(len(sequence))]
    for i in range(len(sequence) - (kmer_size - 1)):
        kmer_count = 0
        start_index = i
        end_index = i + (kmer_size - 1)
        for j in range(i, len(sequence), kmer_size):
            if sequence[i:i+kmer_size] == sequence[j:j+kmer_size]:
                kmer_count += 1
            else:
                break
            end_index = j + (kmer_size)
        for k in range(start_index, min(len(sequence), end_index)):
            max_observed_repeats[k] = max(max_observed_repeats[k], kmer_count)

    return max_observed_repeats

def overlapping_windows(sequence, window_size):
    """
    Implemented from: https://jszym.com/blog/dna_protein_complexity/
    Returns overlapping windows of size `window_size` from sequence `sequence`
    :param sequence: Sequence to extract windows from
    :param window_size: the length of the windows to yield
    """
    for i in range(len(sequence) - window_size + 1):
        yield sequence[i:i+window_size]


def sequence_to_repvec(sequence):
    """
    Implemented from: https://jszym.com/blog/dna_protein_complexity/
    Computes the repetition vector (as seen in Wooton, 1993) from a
    given sequence of a biopolymer with `N` possible residues.

    :param sequence: the nucleotide or protein sequence to generate a repetition vector for.
    """
    encountered_homopolymers = set()
    repvec = []

    for base in sequence:
        if base not in encountered_homopolymers:
            residue_count = sequence.count(base)

            repvec.append(residue_count)

            encountered_homopolymers.add(base)

        if len(encountered_homopolymers) == 4:
            break

    while len(repvec) < 4:
        repvec.append(0)

    return sorted(repvec, reverse=True)


def sequence_entropy(sequence):
    """
    Implemented from: https://jszym.com/blog/dna_protein_complexity/
    Computes the Shannon Entropy of a given sequence. See (Wooton, 1993) for more.

    :param sequence: the nucleotide or protein sequence whose Shannon Entropy is to calculated.
    """
    repvec = sequence_to_repvec(sequence)

    entropy = 0

    for n in repvec:
        # we have to throw away elements whose value is 0 since log(0) is invalid
        if n != 0:
            entropy += -1 * ((n/len(sequence))*math.log((n/len(sequence)), 4))

    return entropy

def small_chunk_stitch(options, file_chunks):
    fasta_handler = PEPPER_VARIANT.FASTA_handler(options.fasta)
    selected_candidate_list_margin = []
    selected_candidate_list_deepvariant = []
    for file_chunk in file_chunks:
        file_name, batch_key = file_chunk
        all_candidates = []
        with h5py.File(file_name, 'r') as hdf5_file:
            if 'predictions' in hdf5_file.keys():
                contigs = hdf5_file['predictions'][batch_key]['contigs'][()]
                positions = hdf5_file['predictions'][batch_key]['positions'][()]
                depths = hdf5_file['predictions'][batch_key]['depths'][()]
                candidates = hdf5_file['predictions'][batch_key]['candidates'][()]
                candidate_frequencies = hdf5_file['predictions'][batch_key]['candidate_frequency'][()]
                base_predictions = hdf5_file['predictions'][batch_key]['base_prediction'][()]
                # type_predictions = hdf5_file['predictions'][batch_key]['type_prediction'][()]

                for i in range(len(contigs)):
                    candidate = str(candidates[i]).strip('][')
                    candidate = candidate.replace(',', ' ')
                    candidate = candidate.split()
                    candidate = [x.strip("'") for x in candidate]
                    candidate_frequency = str(candidate_frequencies[i]).strip('][')
                    candidate_frequency = candidate_frequency.replace(',', ' ')
                    candidate_frequency = candidate_frequency.split()
                    candidate_frequency = [int(x.strip("'")) for x in candidate_frequency]
                    candidate = PEPPER_VARIANT.CandidateImagePrediction(contigs[i].decode('UTF-8'),
                                                                        positions[i],
                                                                        depths[i],
                                                                        candidate,
                                                                        candidate_frequency,
                                                                        base_predictions[i],
                                                                        [])
                    all_candidates.append(candidate)

        for candidate in all_candidates:

            reference_base = fasta_handler.get_reference_sequence(candidate.contig, candidate.position, candidate.position+1).upper()
            reference_upstream = fasta_handler.get_reference_sequence(candidate.contig, candidate.position, candidate.position + 10).upper()
            reference_downstream = fasta_handler.get_reference_sequence(candidate.contig, max(0, candidate.position - 10), candidate.position).upper()

            full_sequence = reference_downstream + reference_upstream
            full_sequence = full_sequence.upper()

            homopolymer_repeats = repeat_annotation(full_sequence, 1)
            # dimer_repeats = repeat_annotation(full_sequence, 2)
            # trimer_repeats = repeat_annotation(full_sequence, 3)

            position_index = len(reference_downstream)
            upward_lookup_index = min(len(homopolymer_repeats), position_index + 4)
            downward_lookup_index = max(0, position_index - 5)
            # max_dimer_count = max(dimer_repeats[downward_lookup_index:upward_lookup_index])
            # max_trimer_count = max(trimer_repeats[downward_lookup_index:upward_lookup_index])
            max_homopolymer_count = max(homopolymer_repeats[downward_lookup_index:upward_lookup_index])

            # upward_entropy = sequence_entropy(fasta_handler.get_reference_sequence(candidate.contig, candidate.position, candidate.position + 10).upper())
            # downward_entropy = sequence_entropy(fasta_handler.get_reference_sequence(candidate.contig, max(0, candidate.position - 10), candidate.position).upper())

            candidate_in_repeat = False
            if max_homopolymer_count >= 5:
                candidate_in_repeat = True

            if reference_base not in ['A', 'C', 'G', 'T']:
                continue

            predicted_genotype = np.argmax(candidate.prediction_base)

            if predicted_genotype == 0:
                genotype = [0, 0]
            elif predicted_genotype == 1:
                genotype = [0, 1]
            else:
                genotype = [1, 1]

            prediction_value = candidate.prediction_base[predicted_genotype]

            # this is for Margin. Only pick SNPs.
            alt_alleles = []
            variant_allele_support = []
            for alt_allele, allele_frequency in zip(candidate.candidates, candidate.candidate_frequency):
                alt_type = alt_allele[0]
                allele = alt_allele[1:]

                allele_list = list(allele)
                valid_allele = True
                for base in allele_list:
                    if base not in ['A', 'C', 'G', 'T']:
                        valid_allele = False

                if not valid_allele:
                    continue
                # only process SNPs for margin
                if alt_type == '1':
                    if predicted_genotype != 0:
                        alt_alleles.append(allele)
                        variant_allele_support.append(allele_frequency)

            if len(alt_alleles) > 0:
                # print(candidate.contig, candidate.position, candidate.position + 1, reference_base, alt_alleles, genotype, candidate.depth, variant_allele_support)
                selected_candidate_list_margin.append((candidate.contig, candidate.position, candidate.position + 1, reference_base, alt_alleles, genotype, candidate.depth, variant_allele_support, prediction_value, candidate.prediction_base))

            ################################################
            # this is candidate finding for DeepVariant
            ################################################
            alt_alleles = []
            variant_allele_support = []
            max_delete_length = 0
            reference_allele = reference_base
            non_alt_predictions = []
            for alt_allele, allele_frequency in zip(candidate.candidates, candidate.candidate_frequency):
                # print("GENERAL: ", candidate.contig, candidate.position, reference_allele, ''.join(alt_allele), candidate.depth, allele_frequency)
                alt_type = alt_allele[0]
                allele = alt_allele[1:]

                allele_list = list(allele)
                valid_allele = True
                for base in allele_list:
                    if base not in ['A', 'C', 'G', 'T']:
                        valid_allele = False

                if not valid_allele:
                    continue

                vaf = float(allele_frequency) / float(candidate.depth)
                non_alt_prediction = max(candidate.prediction_base[1], candidate.prediction_base[2])
                non_alt_predictions.append(non_alt_prediction)
                if alt_type == '1':
                    if not candidate_in_repeat and non_alt_prediction >= options.snp_p_value:
                        # add them to list
                        alt_alleles.append(''.join(alt_allele[1:]))
                        variant_allele_support.append(allele_frequency)
                    # repeat variants
                    elif candidate_in_repeat and non_alt_prediction >= options.snp_p_value_in_lc:
                        alt_alleles.append(''.join(alt_allele[1:]))
                        variant_allele_support.append(allele_frequency)
                    elif 0 < options.report_snp_above_freq <= vaf:
                        alt_alleles.append(''.join(alt_allele[1:]))
                        variant_allele_support.append(allele_frequency)
                elif alt_type == '2':
                    if not candidate_in_repeat and non_alt_prediction >= options.insert_p_value:
                        # add them to list
                        alt_alleles.append(''.join(alt_allele[1:]))
                        variant_allele_support.append(allele_frequency)
                    # repeat variants
                    elif candidate_in_repeat and non_alt_prediction >= options.insert_p_value_in_lc:
                        alt_alleles.append(''.join(alt_allele[1:]))
                        variant_allele_support.append(allele_frequency)
                    elif 0 < options.report_indel_above_freq <= vaf:
                        alt_alleles.append(''.join(alt_allele[1:]))
                        variant_allele_support.append(allele_frequency)
                elif alt_type == '3':
                    if not candidate_in_repeat and non_alt_prediction >= options.delete_p_value:
                        # add them to list
                        alt_alleles.append(reference_allele)
                        reference_allele = ''.join(alt_allele[1:])
                        variant_allele_support.append(allele_frequency)
                    # repeat variants
                    elif candidate_in_repeat and non_alt_prediction >= options.delete_p_value_in_lc:
                        # add them to list
                        alt_alleles.append(reference_allele)
                        reference_allele = ''.join(alt_allele[1:])
                        variant_allele_support.append(allele_frequency)
                    elif 0 < options.report_indel_above_freq <= vaf:
                        alt_alleles.append(''.join(alt_allele[1:]))
                        variant_allele_support.append(allele_frequency)

                predicted_genotype = np.argmax(candidate.prediction_base)
                # print(candidate.contig, candidate.position, reference_base, alt_allele, allele_frequency, candidate.depth, vaf, "PRED: ", non_alt_prediction, "QUAL", non_alt_phred, candidate.prediction_base, candidate_in_repeat)
            if len(alt_alleles) > 0:
                # print("SELECTED", candidate.contig, candidate.position, candidate.position + 1, reference_allele, alt_alleles, genotype, candidate.depth, variant_allele_support, prediction_value, non_alt_phred, candidate_in_repeat)
                selected_candidate_list_deepvariant.append((candidate.contig, candidate.position, candidate.position + len(reference_allele), reference_allele, alt_alleles, genotype, candidate.depth, variant_allele_support, prediction_value, candidate.prediction_base, non_alt_predictions, candidate_in_repeat))
            # print("----------------------------------------------------------------------")

    return selected_candidate_list_margin, selected_candidate_list_deepvariant


def find_candidates(options, input_dir, all_prediction_pair):

    all_selected_candidates_phasing = list()
    all_selected_candidates_variant_calling = list()
    # generate the dictionary in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=options.threads) as executor:
        file_chunks = chunks(all_prediction_pair, max(2, int(len(all_prediction_pair) / options.threads) + 1))
        futures = [executor.submit(small_chunk_stitch, options, file_chunk) for file_chunk in file_chunks]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                positional_candidates_phasing, positional_candidates_variant_calling = fut.result()
                all_selected_candidates_phasing.extend(positional_candidates_phasing)
                all_selected_candidates_variant_calling.extend(positional_candidates_variant_calling)
            else:
                sys.stderr.write("ERROR IN THREAD: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    all_selected_candidates_phasing = sorted(all_selected_candidates_phasing, key=lambda x: (x[0], x[1]))
    all_selected_candidates_variant_calling = sorted(all_selected_candidates_variant_calling, key=lambda x: (x[0], x[1]))

    all_selected_candidates_phasing_positional_dictionary = defaultdict(list)
    all_selected_candidates_phasing_positional_dictionary_alts = defaultdict(list)
    all_selected_candidates_variant_calling_positional_dictionary = defaultdict(list)
    all_selected_candidates_variant_calling_positional_dictionary_alts = defaultdict(list)

    for candidate in all_selected_candidates_phasing:
        ref = candidate[3]
        alt = candidate[4][0]
        if (ref, alt) in all_selected_candidates_phasing_positional_dictionary_alts[(candidate[0], candidate[1])]:
            continue
        all_selected_candidates_phasing_positional_dictionary_alts[(candidate[0], candidate[1])].append((ref, alt))
        all_selected_candidates_phasing_positional_dictionary[(candidate[0], candidate[1])].append(candidate)

    contigs = list()
    for candidate in all_selected_candidates_variant_calling:
        if candidate[0] not in contigs:
            contigs.append(candidate[0])
        ref = candidate[3]
        alt = candidate[4][0]
        if (ref, alt) in all_selected_candidates_variant_calling_positional_dictionary_alts[(candidate[0], candidate[1])]:
            continue
        all_selected_candidates_variant_calling_positional_dictionary_alts[(candidate[0], candidate[1])].append((ref, alt))
        all_selected_candidates_variant_calling_positional_dictionary[(candidate[0], candidate[1])].append(candidate)

    # print(all_selected_candidates_variant_calling_positional_dictionary)
    # print(all_selected_candidates_phasing_positional_dictionary)
    # exit(0)

    # print("SORTED")
    return contigs, all_selected_candidates_phasing_positional_dictionary, all_selected_candidates_variant_calling_positional_dictionary
