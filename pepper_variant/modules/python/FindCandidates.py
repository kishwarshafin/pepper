import h5py
from datetime import datetime
import sys
from os.path import isfile, join
from os import listdir
import re
import time
import pickle
from pepper_variant.build import PEPPER_VARIANT
from pepper_variant.modules.python.CandidateFinder import find_candidates
from pepper_variant.modules.python.VcfWriter import VCFWriter
from pepper_variant.modules.python.ImageGenerationUI import ImageGenerationUtils
from pepper_variant.modules.python.Options import CandidateFinderOptions


def candidates_to_variants(candidates, contig):
    max_h1_prob = 0.0
    max_h2_prob = 0.0
    h1_indx = -1
    h2_indx = -1
    min_pos_start = -1
    max_pos_end = -1
    ref_sequence = ""
    overall_non_ref_prob = -1.0

    for i, candidate in enumerate(candidates):
        pos_start, pos_end, ref, alt, alt_type, depth, read_support, \
            read_support_h0, read_support_h1, read_support_h2, alt_prob_h1, alt_prob_h2, non_ref_prob = candidate

        if overall_non_ref_prob < 0:
            overall_non_ref_prob = non_ref_prob

        overall_non_ref_prob = min(non_ref_prob, overall_non_ref_prob)

        if min_pos_start == -1:
            min_pos_start = pos_start
        if max_pos_end == -1:
            max_pos_end = pos_end

        min_pos_start = min(min_pos_start, pos_start)
        max_pos_end = max(max_pos_end, pos_end)

        if max_pos_end == pos_end:
            ref_sequence = ref

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
    # print(candidates)
    # print(h1_indx, h2_indx)

    selected_alts = []
    selected_dps = []
    selected_gts = []
    selected_ads = []

    other_alts = []
    other_dps = []
    other_gts = []
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
            selected_gts.append(max(alt_prob_h1, alt_prob_h2))
        else:
            other_alts.append(alt)
            other_dps.append(depth)
            other_ads.append(read_support)
            other_gts.append(max(alt_prob_h1, alt_prob_h2))

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
    gts = selected_gts + other_gts
    ads = selected_ads + other_ads
    return contig, min_pos_start, max_pos_end, ref_sequence, alleles, genotype, dps, gts, ads, overall_non_ref_prob


def simplify_variants(variant):
    contig, ref_start, ref_end, ref_seq, alleles, genotype = variant
    if len(alleles) > 1:
        print("ERROR: OBSERVED MORE THAN 1 CANDIDATES AT SITE: ", contig, ref_start, alleles)
        exit(1)
    allele = alleles[0]

    if len(allele) == 1 or len(ref_seq) == 1:
        return [(contig, ref_start, ref_end, ref_seq, alleles, genotype)]

    window_move = min(len(ref_seq), len(allele))
    simplified_variants = []
    for pos in range(ref_start, ref_start + window_move - 1):
        indx = pos - ref_start
        ref_base = ref_seq[indx]
        alt_base = allele[indx]
        if ref_base == alt_base:
            continue
        simplified_variants.append((contig, pos, pos+1, ref_base, [alt_base], genotype))

    ref_out = ref_seq[window_move-1:]
    alt_out = allele[window_move-1:]
    if ref_out != alt_out:
        simplified_variants.append((contig, ref_start+window_move-1, ref_end, ref_seq[window_move-1:], [allele[window_move-1:]], genotype))
    return simplified_variants


def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-3:] == 'hdf']
    return file_paths


def candidate_finder(options, input_dir, output_path):
    all_prediction_files = get_file_paths_from_directory(input_dir)

    all_prediction_pair = []

    for prediction_file in all_prediction_files:
        with h5py.File(prediction_file, 'r') as hdf5_file:
            if 'predictions' in hdf5_file.keys():
                batches = list(hdf5_file['predictions'].keys())
                for batch in batches:
                    all_prediction_pair.append((prediction_file, batch))

    vcf_file_name_full = "PEPPER_VARIANT_FULL"
    vcf_file_name_variant_calling = "PEPPER_VARIANT_OUTPUT_VARIANT_CALLING"
    vcf_file_name_pepper = "PEPPER_VARIANT_OUTPUT_PEPPER"

    local_start_time = time.time()
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: STARTING CANDIDATE FINDING." + "\n")

    contigs, selected_candidates_phasing, selected_candidates_variant_calling = find_candidates(options, input_dir, all_prediction_pair)
    end_time = time.time()

    vcf_file_full = VCFWriter(contigs, options.fasta, options.sample_name, output_path, vcf_file_name_full, vcf_file_name_pepper, vcf_file_name_variant_calling)

    mins = int((end_time - local_start_time) / 60)
    secs = int((end_time - local_start_time)) % 60

    # vcf_file_phasing.write_vcf_records(selected_candidates_phasing, options, calling_mode=0)
    total_variants, total_pepper, total_variant_calling, total_variant_calling_snp, total_variant_calling_indel = vcf_file_full.write_vcf_records(selected_candidates_variant_calling, options)
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING, TOTAL CANDIDATES FOUND: " + str(total_variants) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING, TOTAL VARIANTS IN PEPPER: " + str(total_pepper) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING, TOTAL VARIANTS SELECTED FOR RE-GENOTYPING: " + str(total_variant_calling) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING, TOTAL SNP VARIANTS SELECTED FOR RE-GENOTYPING: " + str(total_variant_calling_snp) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING, TOTAL INDEL VARIANTS SELECTED FOR RE-GENOTYPING: " + str(total_variant_calling_indel) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL TIME SPENT ON CANDIDATE FINDING: " + str(mins) + " Min " + str(secs) + " Sec\n")


def process_candidates(options, input_dir, output_dir):
    output_dir = ImageGenerationUtils.handle_output_directory(output_dir)

    candidate_finder(options,
                     input_dir,
                     output_dir)
