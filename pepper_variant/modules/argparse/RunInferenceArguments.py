def add_run_inference_arguments(parser):
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
        default=512,
        help="Batch size for testing, default is 100. Suggested values: 256/512/1024. Default is 512."
    )
    parser.add_argument(
        "--quantized",
        default=False,
        action='store_true',
        help="PEPPER: Use quantization for inference while on CPU inference mode. Speeds up inference. Default is True."
    )
    parser.add_argument(
        "--no_quantized",
        dest='quantized',
        default=False,
        action='store_false',
        help="Do not use quantization for inference while on CPU inference mode. Speeds up inference."
    )
    parser.add_argument(
        "-g",
        "--gpu",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs for inference. CUDA required. Default is False."
    )
    parser.add_argument(
        "-d_ids",
        "--device_ids",
        type=str,
        required=False,
        default=None,
        help="List of gpu device ids to use for inference. Only used in distributed setting.\n"
             "Example usage: --device_ids 0,1,2 (this will create three callers in id 'cuda:0, cuda:1 and cuda:2'\n"
             "If none then it will use all available devices. Default is None."
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
        "-w",
        "--num_workers",
        type=int,
        required=False,
        default=0,
        help="Number of workers for loading images. Default is 0."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=8,
        help="Total threads. Default is 8."
    )
    parser.add_argument(
        "-hp",
        "--use_hp_info",
        default=False,
        action='store_true',
        help="If set then haplotype-aware mode will be enabled."
    )
    parser.add_argument(
        "--dry",
        default=False,
        action='store_true',
        help="If set then predictions will be imitated from the training set. [DEBUG MODE]"
    )
    profile_group = parser.add_mutually_exclusive_group(required=True)
    profile_group.add_argument("--ont_r9_guppy5_sup",
                               default=False,
                               action='store_true',
                               help="Set to call variants on R9.4.1 Guppy 5+ sup Oxford Nanopore reads.")
    profile_group.add_argument("--ont_r9_guppy4_hac",
                               default=False,
                               action='store_true',
                               help="Set to call variants on R9.4.1 Guppy HAC caller Oxford Nanopore reads.")
    profile_group.add_argument("--ont_r10_q20",
                               default=False,
                               action='store_true',
                               help="Set to call variants on R10.4 Q20 Oxford Nanopore reads.")
    profile_group.add_argument("--hifi",
                               default=False,
                               action='store_true',
                               help="Set to call variants on PacBio HiFi reads.")
    profile_group.add_argument("--clr",
                               default=False,
                               action='store_true',
                               help="Set to call variants on PacBio CLR reads.")
    return parser