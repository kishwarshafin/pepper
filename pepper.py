import argparse
import sys
import torch
from modules.python.TextColor import TextColor
from modules.python.make_images import make_images
from modules.python.call_consensus import call_consensus
from modules.python.perform_stitch import perform_stitch


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s: string holding boolean value
    :return:
    """
    if s.lower() not in {'false', 'true', '1', 't', '0', 'f'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true' or s.lower() == 't' or s.lower() == '1'


def add_make_images_arguments(parser):
    """
    Add arguments to a parser for sub-command "make_images"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        help="BAM file containing mapping between reads and the draft assembly."
    )
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=True,
        help="FASTA file containing the draft assembly."
    )
    parser.add_argument(
        "-r",
        "--region",
        type=str,
        help="Region in [contig_name:start-end] format"
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        default="make_image_output/",
        help="Path to output directory, if it doesn't exist it will be created."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=5,
        help="Number of threads to use. Default is 5."
    )
    return parser


def add_call_consensus_arguments(parser):
    """
    Add arguments to a parser for sub-command "call_consensus"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-i",
        "--image_dir",
        type=str,
        required=True,
        help="Path to directory containing all HDF5 images."
    )
    parser.add_argument(
        "-m",
        "--model_path",
        type=str,
        required=True,
        help="Path to a trained model."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        default='output',
        help="Path to the output directory."
    )
    parser.add_argument(
        "-b",
        "--batch_size",
        type=int,
        required=False,
        default=128,
        help="Batch size for testing, default is 100. Suggested values: 256/512/1024."
    )

    parser.add_argument(
        "-g",
        "--gpu",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs for inference. CUDA required."
    )
    parser.add_argument(
        "-d",
        "--distributed",
        default=True,
        action='store_true',
        help="Distributed gpu inference. This mode will enable one caller per GPU."
    )
    parser.add_argument(
        "-d_ids",
        "--device_ids",
        type=str,
        required=False,
        default=None,
        help="List of gpu device ids to use for inference. Only used in distributed setting.\n"
             "Example usage: --device_ids 0,1,2 (this will create three callers in id 'cuda:0, cuda:1 and cuda:2'\n"
             "If none then it will use all available devices."
    )
    parser.add_argument(
        "-onnx_off",
        "--onnx_off",
        default=True,
        action='store_false',
        help="Turn off cpu acceleration mode (Disabled when GPU is in use)."
    )
    parser.add_argument(
        "-w",
        "--num_workers",
        type=int,
        required=False,
        default=4,
        help="Number of workers for loading images. Default is 4."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=8,
        help="Total threads to be used per caller. A sane value would be num_callers * threads <= total_threads."
    )
    parser.add_argument(
        "-c",
        "--callers",
        type=int,
        required=False,
        default=8,
        help="Total number of callers to spawn while doing CPU inference in distributed mode."
    )
    return parser


def add_stitch_arguments(parser):
    """
    Add arguments to a parser for sub-command "stitch"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-i",
        "--input_dir",
        type=str,
        required=True,
        help="Input dir containing hdf prediction file."
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="Path to output file with an expected prefix (i.e. -o ./outputs/polished_genome)"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=5,
        help="Number of threads."
    )
    return parser


def main():
    """
    Main interface for PEPPER. The submodules supported as of now are these:
    1) Make images
    2) Call consensus
    3) Stitch
    """
    parser = argparse.ArgumentParser(description="PEPPER is a RNN based polisher for polishing ONT-based assemblies. "
                                                 "It works in three steps:\n"
                                                 "1) make_images: This module takes alignment file and coverts them"
                                                 "to HDF5 files containing summary statistics.\n"
                                                 "2) call_consensus: This module takes the summary images and a"
                                                 "trained neural network and generates predictions per base.\n"
                                                 "3) stitch: This module takes the inference files as input and "
                                                 "stitches them to generate a polished assembly.\n",
                                     formatter_class=argparse.RawTextHelpFormatter)

    subparsers = parser.add_subparsers(dest='sub_command')
    subparsers.required = True

    parser_make_images = subparsers.add_parser('make_images', help="Generate images that encode summary statistics "
                                                                   "of reads aligned to an assembly.")
    add_make_images_arguments(parser_make_images)

    parser_call_consensus = subparsers.add_parser('call_consensus', help="Perform inference on generated images using "
                                                                         "a trained model.")
    add_call_consensus_arguments(parser_call_consensus)

    parser_stitch = subparsers.add_parser('stitch', help="Stitch the polished genome to generate a contiguous polished"
                                                         "assembly.")
    add_stitch_arguments(parser_stitch)
    parser_torch_stat = subparsers.add_parser('torch_stat', help="See PyTorch configuration.")

    FLAGS, unparsed = parser.parse_known_args()

    if FLAGS.sub_command == 'make_images':
        sys.stderr.write(TextColor.GREEN + "INFO: MAKE IMAGE MODULE SELECTED\n" + TextColor.END)
        make_images(FLAGS.bam,
                    FLAGS.fasta,
                    FLAGS.region,
                    FLAGS.output_dir,
                    FLAGS.threads)

    elif FLAGS.sub_command == 'call_consensus':
        sys.stderr.write(TextColor.GREEN + "INFO: CALL CONSENSUS MODULE SELECTED\n" + TextColor.END)
        call_consensus(FLAGS.image_dir,
                       FLAGS.model_path,
                       FLAGS.batch_size,
                       FLAGS.num_workers,
                       FLAGS.output_dir,
                       FLAGS.device_ids,
                       FLAGS.gpu,
                       FLAGS.distributed,
                       FLAGS.callers,
                       FLAGS.threads)

    elif FLAGS.sub_command == 'stitch':
        sys.stderr.write(TextColor.GREEN + "INFO: STITCH MODULE SELECTED\n" + TextColor.END)
        perform_stitch(FLAGS.input_dir,
                       FLAGS.output_file,
                       FLAGS.threads)

    elif FLAGS.sub_command == 'torch_stat':
        sys.stderr.write(TextColor.YELLOW + "TORCH VERSION: " + TextColor.END + str(torch.__version__) + "\n\n")
        sys.stderr.write(TextColor.YELLOW + "PARALLEL CONFIG:\n" + TextColor.END)
        print(torch.__config__.parallel_info())
        sys.stderr.write(TextColor.YELLOW + "BUILD CONFIG:\n" + TextColor.END)
        print(*torch.__config__.show().split("\n"), sep="\n")

        sys.stderr.write(TextColor.GREEN + "CUDA AVAILABLE: " + TextColor.END + str(torch.cuda.is_available()) + "\n")
        sys.stderr.write(TextColor.GREEN + "GPU DEVICES: " + TextColor.END + str(torch.cuda.device_count()) + "\n")


if __name__ == '__main__':
    main()
