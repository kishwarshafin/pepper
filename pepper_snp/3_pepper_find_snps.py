import h5py
import argparse
import sys
from modules.python.TextColor import TextColor
from modules.python.StitchV2 import create_consensus_sequence
from modules.python.VcfWriter import VCFWriter
from os.path import isfile, join
from os import listdir


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-3:] == 'hdf']
    return file_paths


def perform_stitch(hdf_file_path, reference_path, output_path, threads, sample_name, p_threshold):
    all_prediction_files = get_file_paths_from_directory(hdf_file_path)

    all_contigs = set()
    for prediction_file in all_prediction_files:
        with h5py.File(prediction_file, 'r') as hdf5_file:
            contigs = list(hdf5_file['predictions'].keys())
            all_contigs.update(contigs)
    all_contigs = list(all_contigs)

    vcf_file = VCFWriter(reference_path, sample_name, output_path, all_contigs)

    for contig in all_contigs:
        sys.stderr.write(TextColor.YELLOW + "INFO: PROCESSING CONTIG: " + contig + "\n" + TextColor.END)

        all_chunk_keys = list()
        for prediction_file in all_prediction_files:
            with h5py.File(prediction_file, 'r') as hdf5_file:
                chunk_keys = sorted(hdf5_file['predictions'][contig].keys())
                for chunk_key in chunk_keys:
                    all_chunk_keys.append((prediction_file, chunk_key))

        all_candidates, reference_dict, positions = create_consensus_sequence(contig,
                                                                              all_chunk_keys,
                                                                              threads,
                                                                              p_threshold)
        sys.stderr.write(TextColor.BLUE + "INFO: FINISHED PROCESSING " + contig + ", TOTAL CANDIDATES FOUND: "
                         + str(len(all_candidates)) + ".\n" + TextColor.END)
        vcf_file.write_vcf_records(contig, all_candidates, reference_dict, positions)


    hdf5_file.close()


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser(description="3_pepper_stitch.py performs the final stitching to generate  "
                                                 "the polished sequences.")
    parser.add_argument(
        "-i",
        "--input_hdf",
        type=str,
        required=True,
        help="Input hdf prediction file."
    )
    parser.add_argument(
        "-r",
        "--input_reference",
        type=str,
        required=True,
        help="Input reference/assembly file."
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        type=str,
        required=True,
        help="Name of the sample."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory."
    )
    parser.add_argument(
        "-p",
        "--probability_threshold",
        type=float,
        required=False,
        default=0.3,
        help="Threshold value for reporting SNPs."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=5,
        help="Number of threads."
    )

    FLAGS, unparsed = parser.parse_known_args()
    perform_stitch(FLAGS.input_hdf,
                   FLAGS.input_reference,
                   FLAGS.output_dir,
                   FLAGS.threads,
                   FLAGS.sample_name,
                   FLAGS.probability_threshold)
