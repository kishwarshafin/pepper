import argparse
import sys
import torch
from pepper_variant.modules.argparse.SetParameters import set_parameters
from pepper_variant.modules.python.ImageGenerationUI import ImageGenerationUtils
from pepper_variant.modules.python.TrainModule import train_pepper_model
from pepper_variant.modules.python.TestModule import do_test
from datetime import datetime
from pepper_variant.modules.argparse.MakeImagesArguments import add_make_images_arguments
from pepper.version import __version__


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s: string holding boolean value
    :return:
    """
    if s.lower() not in {'false', 'true', '1', 't', '0', 'f'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true' or s.lower() == 't' or s.lower() == '1'


def add_make_train_images_arguments(parser):
    """
    Add arguments to a parser for sub-command "make_images"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-tv",
        "--truth_vcf",
        type=str,
        required=True,
        default=None,
        help="VCF file containing the truth variants [Usually from GIAB]."
    )
    parser.add_argument(
        "-p",
        "--random_draw_probability",
        required=True,
        type=float,
        help="Random draw probabilty will set the percentage of ref examples to keep in the training set. (Recommended: 0.33)"
    )
    return parser


def add_train_model_arguments(parser):
    parser.add_argument(
        "-train",
        "--train_image_dir",
        type=str,
        required=True,
        help="Training data directory containing HDF files."
    )
    parser.add_argument(
        "-test",
        "--test_image_dir",
        type=str,
        required=True,
        help="Testing data directory containing HDF files."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=False,
        default='./pepper_model_out',
        help="Path to output directory."
    )
    parser.add_argument(
        "-bs",
        "--batch_size",
        type=int,
        required=False,
        default=64,
        help="Batch size for training, default is 64."
    )

    parser.add_argument(
        "--test_batch_size",
        type=int,
        required=False,
        default=64,
        help="Batch size for test, default is 64."
    )
    parser.add_argument(
        "-s",
        "--step_size",
        type=int,
        required=False,
        default=100,
        help="Step size that counts toward each epoch."
    )
    parser.add_argument(
        "-e",
        "--epoch_size",
        type=int,
        required=False,
        default=10,
        help="Epoch size for training iteration."
    )
    parser.add_argument(
        "-w",
        "--num_workers",
        type=int,
        required=False,
        default=16,
        help="Number of data loaders to use."
    )
    parser.add_argument(
        "-lr",
        "--learning_rate",
        type=float,
        required=False,
        default=0.0001,
        help="Learning rate. Default: 0.0001."
    )
    parser.add_argument(
        "-wd",
        "--weight_decay",
        type=float,
        required=False,
        default=0.000001,
        help="Weight decay. Default: 0.000001."
    )
    parser.add_argument(
        "-g",
        "--gpu",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs for inference. CUDA required."
    )
    parser.add_argument(
        "-per_gpu",
        "--callers_per_gpu",
        type=int,
        required=False,
        default=4,
        help="Number of callers to initialize per GPU, on a 11GB GPU, you can go up to 10. Default is 4."
    )
    parser.add_argument(
        "-dx",
        "--distributed_off",
        default=False,
        action='store_true',
        help="Turn off distributed inference. This mode will disable the use of multiple callers."
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
        "-rm",
        "--retrain_model",
        default=False,
        action='store_true',
        help="If set then retrain a pre-trained mode."
    )
    parser.add_argument(
        "-rmp",
        "--retrain_model_path",
        type=str,
        default=False,
        help="Path to the model that will be retrained."
    )
    parser.add_argument(
        "-hp",
        "--use_hp_info",
        default=False,
        action='store_true',
        help="If set then haplotype-aware mode will be enabled."
    )

    return parser


def add_test_model_arguments(parser):
    parser.add_argument(
        "--test_image_dir",
        type=str,
        required=True,
        help="Training data description csv file."
    )
    parser.add_argument(
        "--model_path",
        type=str,
        required=True,
        default='./model',
        help="Path of the model to load and retrain"
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        required=False,
        default=100,
        help="Batch size for training, default is 100."
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        required=False,
        default=40,
        help="Epoch size for training iteration."
    )
    parser.add_argument(
        "--gpu",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs. CUDA required."
    )
    parser.add_argument(
        "--print_details",
        action='store_true',
        default=False,
        help="If true then prints debug messages."
    )
    parser.add_argument(
        "-hp",
        "--use_hp_info",
        default=False,
        action='store_true',
        help="If set then haplotype-aware mode will be enabled."
    )
    return parser


def add_run_hyperband_arguments(parser):
    parser.add_argument(
        "--train_image_dir",
        type=str,
        required=True,
        help="Training data directory containing HDF files."
    )
    parser.add_argument(
        "--test_image_dir",
        type=str,
        required=True,
        help="Training data directory containing HDF files."
    )
    parser.add_argument(
        "--gpu",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs. CUDA required."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        default='./hyperband_output/',
        help="Directory to save the model"
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        required=False,
        default=100,
        help="Batch size for training, default is 100."
    )
    parser.add_argument(
        "--max_epochs",
        type=int,
        required=False,
        default=5,
        help="Max epoch size for training iteration."
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        required=False,
        default=40,
        help="Number of data loader workers."
    )
    return parser


def main():
    """
    Main interface for PEPPER. The submodules supported as of now are these:
    1) Make images
    2) Call consensus
    3) Stitch
    """
    parser = argparse.ArgumentParser(description="PEPPER is a RNN based polisher for polishing ONT-based assemblies.\n"
                                                 "PEPPER_train provides an interface to train a PEPPER model.\n"
                                                 "Available sub-modules:\n"
                                                 "1) make_train_images: Generate training samples.\n"
                                                 "2) train_model: Train a model using train and test files\n"
                                                 "3) test_model: Test a model on a set of training samples.\n"
                                                 "4) run_hyperband: Hyper-parameter tuning [in development].\n",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )

    subparsers = parser.add_subparsers(dest='sub_command')

    parser_make_images = subparsers.add_parser('make_train_images', help="Generate training samples")
    add_make_images_arguments(parser_make_images)
    add_make_train_images_arguments(parser_make_images)

    parser_train_model = subparsers.add_parser('train_model', help="Train a model.")
    add_train_model_arguments(parser_train_model)

    parser_test_model = subparsers.add_parser('test_model', help="Test a pre-trained model.")
    add_test_model_arguments(parser_test_model)

    # parser_test_model = subparsers.add_parser('run_hyperband', help="Run hyperband to find best set of parameters.")
    # add_test_model_arguments(parser_test_model)

    parser_torch_stat = subparsers.add_parser('torch_stat', help="See PyTorch configuration.")

    options, unparsed = parser.parse_known_args()
    options.train_mode = True

    if options.sub_command == 'make_train_images':
        options = set_parameters(options)
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MAKE TRAIN IMAGE MODULE SELECTED\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: DOWNSAMPLE RATE:\t" + str(options.downsample_rate) + "\n")
        options.image_output_directory = options.output_dir
        if not options.region_bed:
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: HIGH CONFIDENCE REGIONS NOT PROVIDED. PLEASE PROVIDE A BED FILE USING -rb.\n")
        ImageGenerationUtils.generate_images(options)
    elif options.sub_command == 'train_model':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TRAIN MODEL MODULE SELECTED\n")
        train_pepper_model(options)

    elif options.sub_command == 'test_model':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TEST MODEL MODULE SELECTED\n")
        do_test(options.test_image_dir,
                options.use_hp_info,
                options.batch_size,
                options.gpu,
                options.num_workers,
                options.model_path)

    elif options.sub_command == 'torch_stat':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TORCH VERSION: " + str(torch.__version__) + "\n\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PARALLEL CONFIG:\n")
        print(torch.__config__.parallel_info())
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: BUILD CONFIG:\n")
        print(*torch.__config__.show().split("\n"), sep="\n")

        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CUDA AVAILABLE: " + str(torch.cuda.is_available()) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: GPU DEVICES: " + str(torch.cuda.device_count()) + "\n")

    elif options.version is True:
        print("PEPPER VERSION: ", __version__)

    else:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: NO SUBCOMMAND SELECTED. PLEASE SELECT ONE OF THE AVAIABLE SUB-COMMANDS.\n")
        parser.print_help()


if __name__ == '__main__':
    main()
