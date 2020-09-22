import h5py
import sys
import re
import time
from pepper.modules.python.Stitch import create_consensus_sequence
from os.path import isfile, join
from pathlib import Path
from os import listdir
from datetime import datetime


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


def number_key(name):
    """
    Sorting: https://stackoverflow.com/questions/4287209/sort-list-of-strings-by-integer-suffix-in-python
    :param name:
    :return:
    """
    parts = re.findall('[^0-9]+|[0-9]+', name)
    L = []
    for part in parts:
        try:
            L.append(int(part))
        except ValueError:
            L.append(part)
    return L


def perform_stitch(hdf_file_path, output_path, threads):
    all_prediction_files = get_file_paths_from_directory(hdf_file_path)
    all_contigs = set()
    time.sleep(5)
    for prediction_file in sorted(all_prediction_files):
        with h5py.File(prediction_file, 'r') as hdf5_file:
            contigs = list(hdf5_file['predictions'].keys())
            all_contigs.update(contigs)

    output_path = output_path + '_pepper_polished.fa'
    output_directory = Path(output_path).resolve().parents[0]
    output_directory.mkdir(parents=True, exist_ok=True)

    consensus_fasta_file = open(output_path, 'w')

    for contig in sorted(all_contigs, key=natural_key):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PROCESSING CONTIG: " + contig + "\n")
        sys.stderr.flush()
        # get all the chunk keys from all the files
        all_chunk_keys = list()
        for prediction_file in all_prediction_files:
            with h5py.File(prediction_file, 'r') as hdf5_file:
                if contig not in hdf5_file['predictions'].keys():
                    continue
                chunk_keys = sorted(hdf5_file['predictions'][contig].keys())
                for chunk_key in chunk_keys:
                    contig_start = hdf5_file['predictions'][contig][chunk_key]['contig_start'][()]
                    contig_end = hdf5_file['predictions'][contig][chunk_key]['contig_end'][()]
                    all_chunk_keys.append((prediction_file, chunk_key, contig_start, contig_end))

        consensus_sequence = create_consensus_sequence(contig, all_chunk_keys, threads)

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PROCESSING " + contig + ", POLISHED SEQUENCE LENGTH: "
                         + str(len(consensus_sequence)) + ".\n")

        # TODO: I should write a FASTA handler here. This is too sloppy.
        if consensus_sequence is not None and len(consensus_sequence) > 0:
            consensus_fasta_file.write('>' + contig + "\n")
            consensus_fasta_file.write(consensus_sequence+"\n")

    hdf5_file.close()
