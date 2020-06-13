import h5py
import sys
from datetime import datetime
import time
from pepper_snp.modules.python.CandidateFinder import find_SNP_candidates
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


def find_candidates(input_dir, reference_path, output_dir, threads, sample_name, p_threshold):
    start_time = time.time()
    output_dir = UserInterfaceSupport.handle_output_directory(output_dir)
    all_prediction_files = get_file_paths_from_directory(input_dir)

    all_contigs = set()
    for prediction_file in all_prediction_files:
        with h5py.File(prediction_file, 'r') as hdf5_file:
            if 'predictions' in hdf5_file.keys():
                contigs = list(hdf5_file['predictions'].keys())
                all_contigs.update(contigs)
    all_contigs = list(all_contigs)

    vcf_file = VCFWriter(reference_path, sample_name, output_dir, all_contigs)

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

        all_candidates, reference_dict, positions = find_SNP_candidates(contig,
                                                                        all_chunk_keys,
                                                                        threads,
                                                                        p_threshold)

        vcf_file.write_vcf_records(contig, all_candidates, reference_dict, positions)

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING " + contig + ", TOTAL CANDIDATES FOUND: "
                         + str(len(all_candidates)) + ".\n")
        end_time = time.time()
        mins = int((end_time - start_time) / 60)
        secs = int((end_time - start_time)) % 60
        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] TOTAL ELAPSED TIME FOR CANDIDATE FINDING: " + str(mins) + " Min " + str(secs) + " Sec\n")

    hdf5_file.close()

