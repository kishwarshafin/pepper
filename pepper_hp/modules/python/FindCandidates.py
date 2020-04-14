import h5py
from datetime import datetime
import sys
from os.path import isfile, join
from os import listdir
import re
from pepper_hp.modules.python.CandidateFinder import find_candidates
from pepper_hp.modules.python.VcfWriter import VCFWriter
from pepper_hp.modules.python.ImageGenerationUI import UserInterfaceSupport


def candidates_to_variants(candidates, contig):
    alleles = []
    ref_start = candidates[0][0]
    ref_end = candidates[0][1]
    ref_seq = candidates[0][2]

    if len(candidates) > 1:
        if candidates[1][1] > candidates[0][1]:
            ref_end = candidates[1][1]
            ref_seq = candidates[1][2]

    for candidate in candidates:
        changed_candidate = list(candidate)
        if candidate[1] < ref_end:
            bases_needed = ref_end - candidate[1]
            ref_suffix = ref_seq[-bases_needed:]
            changed_candidate[1] = ref_end
            changed_candidate[2] = ref_seq
            changed_candidate[3] = candidate[3] + ref_suffix
        alleles.append(changed_candidate[3])

    genotype = [0, 0]
    if len(alleles) == 1:
        genotype = [0, 1]
    elif alleles[0] == alleles[1]:
        alleles = list(set(alleles))
        genotype = [1, 1]
    else:
        genotype = [1, 2]

    return contig, ref_start, ref_end, ref_seq, alleles, genotype


def write_vcf(contig, all_candidate_positions, candidate_positional_map, vcf_file):
    # print(candidate_map)
    # candidate_map = {2931716: {(2931716, 2931719, 'CTT', 'C', 1, 'DEL'), (2931716, 2931718, 'CT', 'C', 2, 'DEL')}}
    for pos in sorted(all_candidate_positions):
        candidates = list()
        if pos in candidate_positional_map:
            candidates.append(candidate_positional_map[pos])

        if len(candidates) == 0:
            continue

        if len(candidates) > 2:
            print("ERROR: OBSERVED MORE THAN 2 CANDIDATES AT SITE: ", contig, pos, candidates)
            exit(1)

        variant = candidates_to_variants(list(candidates), contig)
        vcf_file.write_vcf_records(variant)


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


def candidate_finder(input_dir, reference_file, sample_name, output_path, threads):
    all_prediction_files = get_file_paths_from_directory(input_dir)

    all_contigs = set()
    for prediction_file in all_prediction_files:
        with h5py.File(prediction_file, 'r') as hdf5_file:
            if 'predictions' in hdf5_file.keys():
                contigs = list(hdf5_file['predictions'].keys())
                all_contigs.update(contigs)
    all_contigs = sorted(all_contigs, key=natural_key)

    vcf_file = VCFWriter(reference_file, all_contigs, sample_name, output_path, "candidates_as_variants")

    for contig in sorted(all_contigs, key=natural_key):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESSING CONTIG: " + contig + "\n")

        all_chunk_keys = list()
        for prediction_file in all_prediction_files:
            with h5py.File(prediction_file, 'r') as hdf5_file:
                if 'predictions' in hdf5_file.keys():
                    chunk_keys = sorted(hdf5_file['predictions'][contig].keys())
                    for chunk_key in chunk_keys:
                        all_chunk_keys.append((prediction_file, chunk_key))

        all_candidate_positions, candidate_positional_map = \
            find_candidates(input_dir,  reference_file, contig, all_chunk_keys, threads)

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING " + contig + ", TOTAL CANDIDATES FOUND: "
                         + str(len(all_candidate_positions)) + ".\n")

        write_vcf(contig, all_candidate_positions, candidate_positional_map, vcf_file)

    hdf5_file.close()


def process_candidates(input_dir, reference, sample_name, output_dir, threads):
    output_dir = UserInterfaceSupport.handle_output_directory(output_dir)
    candidate_finder(input_dir,
                     reference,
                     sample_name,
                     output_dir,
                     threads)
