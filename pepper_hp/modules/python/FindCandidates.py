import h5py
from datetime import datetime
import sys
from os.path import isfile, join
from os import listdir
import re
import time
from pepper_hp.build import PEPPER_HP
from pepper_hp.modules.python.CandidateFinder import find_candidates, find_candidates_ccs
from pepper_hp.modules.python.VcfWriter import VCFWriter
from pepper_hp.modules.python.ImageGenerationUI import UserInterfaceSupport
from pepper_hp.modules.python.Options import CandidateFinderOptions
from pepper_hp.modules.python.ExcludeContigs import EXCLUDED_HUMAN_CONTIGS


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


def candidate_finder(input_dir, reference_file, bam_file, sample_name, output_path, threads, split_candidates, set_profile):
    hp_tags = [0]

    if split_candidates:
        hp_tags = [1, 2]

    for haplotag in hp_tags:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESSING HAPLOTAG: " + str(haplotag) + "\n")
        all_prediction_files = get_file_paths_from_directory(input_dir)

        all_contigs = set()
        for prediction_file in all_prediction_files:
            with h5py.File(prediction_file, 'r') as hdf5_file:
                if 'predictions' in hdf5_file.keys():
                    contigs = list(hdf5_file['predictions'].keys())
                    all_contigs.update(contigs)
        all_contigs = sorted(all_contigs, key=natural_key)

        vcf_file = VCFWriter(reference_file, all_contigs, sample_name, output_path, "PEPPER_HP_OUTPUT_" + str(haplotag))

        for contig in sorted(all_contigs, key=natural_key):
            if len(hp_tags) > 1:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESSING CONTIG: " + contig + " WITH HP-TAG: " + str(haplotag) + "\n")
            else:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESSING CONTIG: " + contig + "\n")

            local_start_time = time.time()
            all_chunk_keys = list()
            for prediction_file in all_prediction_files:
                with h5py.File(prediction_file, 'r') as hdf5_file:
                    if 'predictions' in hdf5_file.keys():
                        if contig in hdf5_file['predictions'].keys():
                            chunk_keys = sorted(hdf5_file['predictions'][contig].keys())
                            for chunk_key in chunk_keys:
                                all_chunk_keys.append((prediction_file, chunk_key))

            selected_candidates = find_candidates(input_dir, reference_file, bam_file, contig, all_chunk_keys, threads, haplotag, set_profile)

            end_time = time.time()
            mins = int((end_time - local_start_time) / 60)
            secs = int((end_time - local_start_time)) % 60
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING " + contig + ", TOTAL CANDIDATES FOUND: "
                             + str(len(selected_candidates)) + " TOTAL TIME SPENT: " + str(mins) + " Min " + str(secs) + " Sec\n")

            vcf_file.write_vcf_records(selected_candidates)

        hdf5_file.close()


def candidate_finder_ccs(reference_file, bam_file, sample_name, output_path, threads, split_candidates):
    hp_tags = [0]

    if split_candidates:
        hp_tags = [1, 2]

    for haplotag in hp_tags:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESSING HAPLOTAG: " + str(haplotag) + "\n")
        fasta_handler = PEPPER_HP.FASTA_handler(reference_file)
        bam_handler = PEPPER_HP.BAM_handler(bam_file)
        bam_contigs = bam_handler.get_chromosome_sequence_names()
        fasta_contigs = fasta_handler.get_chromosome_names()
        common_contigs = list(set(fasta_contigs) & set(bam_contigs))
        common_contigs = list(set(common_contigs) - set(EXCLUDED_HUMAN_CONTIGS))
        all_contigs = common_contigs

        vcf_file = VCFWriter(reference_file, all_contigs, sample_name, output_path, "PEPPER_HP_OUTPUT_" + str(haplotag))

        for contig in sorted(all_contigs, key=natural_key):
            if len(hp_tags) > 1:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESSING CONTIG: " + contig + " WITH HP-TAG: " + str(haplotag) + "\n")
            else:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESSING CONTIG: " + contig + "\n")
            local_start_time = time.time()
            all_chunk_keys = list()

            interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(contig) - 1)
            interval_size = interval_end - interval_start

            # this is the interval size each of the process is going to get which is 10^6
            # I will split this into 10^4 size inside the worker process
            all_intervals = []
            for pos in range(interval_start, interval_end, 100000):
                pos_start = max(interval_start, pos)
                pos_end = min(interval_end, pos + 100000)

                all_intervals.append((contig, pos_start, pos_end))

            selected_candidates = find_candidates_ccs(all_intervals, reference_file, bam_file, contig, threads, haplotag)

            end_time = time.time()
            mins = int((end_time - local_start_time) / 60)
            secs = int((end_time - local_start_time)) % 60
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING " + contig + ", TOTAL CANDIDATES FOUND: "
                             + str(len(selected_candidates)) + " TOTAL TIME SPENT: " + str(mins) + " Min " + str(secs) + " Sec\n")

            vcf_file.write_vcf_records(selected_candidates)


def process_candidates(input_dir, reference, bam_file, sample_name, output_dir, threads, split_candidates, set_profile):
    output_dir = UserInterfaceSupport.handle_output_directory(output_dir)

    candidate_finder(input_dir,
                     reference,
                     bam_file,
                     sample_name,
                     output_dir,
                     threads,
                     split_candidates,
                     set_profile)


def process_candidates_ccs(reference, bam_file, sample_name, output_dir, threads, split_candidates):
    output_dir = UserInterfaceSupport.handle_output_directory(output_dir)
    candidate_finder_ccs(reference,
                         bam_file,
                         sample_name,
                         output_dir,
                         threads,
                         split_candidates)
