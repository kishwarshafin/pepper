import argparse
import sys
import torch
from pepper.version import __version__
from pepper.modules.python.polish import polish
from datetime import datetime
from pepper.modules.python.make_images import make_images
from pepper.modules.python.call_consensus import call_consensus
from pepper.modules.python.perform_stitch import perform_stitch
from pepper.modules.python.download_model import download_models


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s: string holding boolean value
    :return:
    """
    if s.lower() not in {'false', 'true', '1', 't', '0', 'f'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true' or s.lower() == 't' or s.lower() == '1'


def add_polish_arguments(parser):
    """
    Add arguments to a parser for sub-command "polish"
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
        "-m",
        "--model_path",
        type=str,
        required=True,
        help="Path to a trained model."
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
        help="Number of threads to use. Default is 5."
    )
    parser.add_argument(
        "-r",
        "--region",
        type=str,
        help="Region in [contig_name:start-end] format"
    )
    parser.add_argument(
        "-bs",
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
        "-w",
        "--num_workers",
        type=int,
        required=False,
        default=4,
        help="Number of workers for loading images. Default is 4."
    )
    return parser


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
        "-bs",
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
        help="Number of threads."
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


def add_download_models_arguments(parser):
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to directory where models will be saved."
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
    parser.add_argument(
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )

    subparsers = parser.add_subparsers(dest='sub_command')
    # subparsers.required = True

    parser_polish = subparsers.add_parser('polish', help="Run the polishing pipeline. This will run "
                                                         "make images-> inference -> stitch one after another.\n"
                                                         "The outputs of each step can be run separately using\n"
                                                         "the appropriate sub-command.")
    add_polish_arguments(parser_polish)

    parser_make_images = subparsers.add_parser('make_images', help="Generate images that encode summary statistics "
                                                                   "of reads aligned to an assembly.")
    add_make_images_arguments(parser_make_images)

    parser_call_consensus = subparsers.add_parser('call_consensus', help="Perform inference on generated images using "
                                                                         "a trained model.")
    add_call_consensus_arguments(parser_call_consensus)

    parser_stitch = subparsers.add_parser('stitch', help="Stitch the polished genome to generate a contiguous polished"
                                                         "assembly.")
    add_stitch_arguments(parser_stitch)

    parser_download_model = subparsers.add_parser('download_models', help="Download available models.")
    add_download_models_arguments(parser_download_model)

    parser_torch_stat = subparsers.add_parser('torch_stat', help="See PyTorch configuration.")
    parser_torch_stat = subparsers.add_parser('version', help="Show program version.")

    FLAGS, unparsed = parser.parse_known_args()

    if FLAGS.sub_command == 'polish':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]  INFO: POLISH MODULE SELECTED\n")
        polish(FLAGS.bam,
               FLAGS.fasta,
               FLAGS.output_file,
               FLAGS.threads,
               FLAGS.region,
               FLAGS.model_path,
               FLAGS.batch_size,
               FLAGS.gpu,
               FLAGS.device_ids,
               FLAGS.num_workers)

    elif FLAGS.sub_command == 'make_images':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MAKE IMAGE MODULE SELECTED\n")
        make_images(FLAGS.bam,
                    FLAGS.fasta,
                    FLAGS.region,
                    FLAGS.output_dir,
                    FLAGS.threads)

    elif FLAGS.sub_command == 'call_consensus':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]  INFO: CALL CONSENSUS MODULE SELECTED\n")
        call_consensus(FLAGS.image_dir,
                       FLAGS.model_path,
                       FLAGS.batch_size,
                       FLAGS.num_workers,
                       FLAGS.output_dir,
                       FLAGS.device_ids,
                       FLAGS.gpu,
                       FLAGS.threads)

    elif FLAGS.sub_command == 'stitch':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: STITCH MODULE SELECTED\n")
        perform_stitch(FLAGS.input_dir,
                       FLAGS.output_file,
                       FLAGS.threads)

    elif FLAGS.sub_command == 'download_models':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: DOWNLOAD MODELS SELECTED\n")
        download_models(FLAGS.output_dir)

    elif FLAGS.sub_command == 'torch_stat':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TORCH VERSION: " + str(torch.__version__) + "\n\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PARALLEL CONFIG:\n")
        print(torch.__config__.parallel_info())
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: BUILD CONFIG:\n")
        print(*torch.__config__.show().split("\n"), sep="\n")

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CUDA AVAILABLE: " + str(torch.cuda.is_available()) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: GPU DEVICES: " + str(torch.cuda.device_count()) + "\n")

    elif FLAGS.version is True:
        print("PEPPER VERSION: ", __version__)

    else:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: NO SUBCOMMAND SELECTED. PLEASE SELECT ONE OF THE AVAIABLE SUB-COMMANDS.\n")
        parser.print_help()


if __name__ == '__main__':
    main()
