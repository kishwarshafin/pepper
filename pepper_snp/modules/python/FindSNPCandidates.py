import h5py
import sys
from datetime import datetime
import time
from pepper_snp.modules.python.CandidateFinder import find_candidates
from pepper_snp.modules.python.VcfWriter import VCFWriter
from pepper_snp.modules.python.ImageGenerationUI import UserInterfaceSupport
from os.path import isfile, join
from os import listdir
import re


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


def snp_candidate_finder(input_dir, reference_file, bam_file, sample_name, output_path, threads, set_profile):
    all_prediction_files = get_file_paths_from_directory(input_dir)

    all_contigs = set()
    for prediction_file in all_prediction_files:
        with h5py.File(prediction_file, 'r') as hdf5_file:
            if 'predictions' in hdf5_file.keys():
                contigs = list(hdf5_file['predictions'].keys())
                all_contigs.update(contigs)
    all_contigs = sorted(all_contigs, key=natural_key)

    vcf_file = VCFWriter(reference_file, all_contigs, sample_name, output_path, "PEPPER_SNP_OUTPUT")

    for contig in sorted(all_contigs, key=natural_key):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESSING CONTIG: " + contig + "\n")

        all_chunk_keys = list()
        for prediction_file in all_prediction_files:
            with h5py.File(prediction_file, 'r') as hdf5_file:
                if 'predictions' in hdf5_file.keys():
                    if contig in hdf5_file['predictions'].keys():
                        chunk_keys = sorted(hdf5_file['predictions'][contig].keys())
                        for chunk_key in chunk_keys:
                            all_chunk_keys.append((prediction_file, chunk_key))

        selected_candidates = find_candidates(reference_file, bam_file, contig, all_chunk_keys, threads, set_profile)

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING " + contig + ", TOTAL CANDIDATES FOUND: "
                         + str(len(selected_candidates)) + ".\n")

        vcf_file.write_vcf_records(selected_candidates)

    hdf5_file.close()

