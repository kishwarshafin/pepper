import h5py
import argparse
import sys
from modules.python.TextColor import TextColor
from modules.python.StitchV2 import create_consensus_sequence


def perform_stitch(hdf_file_path, output_path, threads):
    with h5py.File(hdf_file_path, 'r') as hdf5_file:
        contigs = list(hdf5_file['predictions'].keys())

    consensus_fasta_file_h1 = open(output_path+'consensus_h1.fa', 'w')
    consensus_fasta_file_h2 = open(output_path+'consensus_h2.fa', 'w')
    for contig in contigs:
        sys.stderr.write(TextColor.YELLOW + "INFO: PROCESSING CONTIG: " + contig + "\n" + TextColor.END)

        with h5py.File(hdf_file_path, 'r') as hdf5_file:
            chunk_keys = sorted(hdf5_file['predictions'][contig].keys())

        consensus_sequence_h1, consensus_sequence_h2 = create_consensus_sequence(hdf_file_path, contig, chunk_keys, threads)
        sys.stderr.write(TextColor.BLUE + "INFO: FINISHED PROCESSING " + contig + ", POLISHED SEQUENCE LENGTH: "
                         + str(len(consensus_sequence_h1)) + " " + str(len(consensus_sequence_h2))
                         + ".\n" + TextColor.END)

        # TODO: I should write a FASTA handler here. This is too sloppy.
        if consensus_sequence_h1 is not None and len(consensus_sequence_h1) > 0:
            consensus_fasta_file_h1.write('>' + contig + "_h1\n")
            consensus_fasta_file_h1.write(consensus_sequence_h1+"\n")

        if consensus_sequence_h2 is not None and len(consensus_sequence_h2) > 0:
            consensus_fasta_file_h2.write('>' + contig + "_h2\n")
            consensus_fasta_file_h2.write(consensus_sequence_h2+"\n")

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
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=5,
        help="Number of threads."
    )

    FLAGS, unparsed = parser.parse_known_args()
    perform_stitch(FLAGS.input_hdf, FLAGS.output_dir, FLAGS.threads)
