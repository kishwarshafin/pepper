def add_call_variant_arguments(parser):
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
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory."
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        type=str,
        required=True,
        help="Name of the sample."
    )
    parser.add_argument(
        "-t",
        "--threads",
        required=True,
        type=int,
        help="Number of threads to use."
    )
    # Image generation optional parameters
    parser.add_argument(
        "-d",
        "--downsample_rate",
        type=float,
        default=1.0,
        help="Downsample rate of reads while generating images. Default is 1.0"
    )
    parser.add_argument(
        "-r",
        "--region",
        type=str,
        default=None,
        help="Region in [contig_name:start-end] format. Default is None."
    )
    parser.add_argument(
        "--region_size",
        type=int,
        required=False,
        default=100000,
        help="Region size in bp used to chunk the genome. Default is 100000."
    )
    parser.add_argument(
        "--region_bed",
        "-rb",
        type=str,
        required=False,
        default=None,
        help="Bed file to process regions only in the bed region. Default is None."
    )
    parser.add_argument(
        "-hp",
        "--use_hp_info",
        default=False,
        action='store_true',
        help="If true use HP information for variant calling. Default is False."
    )
    parser.add_argument(
        "--include_supplementary",
        default=False,
        action='store_true',
        help="If true then supplementary reads will be used. Default is False."
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        required=False,
        default=5,
        help="Minimum mapping quality for read to be considered valid. Default is 5"
    )
    parser.add_argument(
        "--min_baseq",
        type=int,
        required=False,
        default=1,
        help="Minimum base quality for base to be considered valid. Default is 1"
    )
    parser.add_argument(
        "--snp_frequency",
        type=float,
        required=False,
        default=0.10,
        help="Minimum SNP frequency for a site to be considered to have a variant. Default is 0.10"
    )
    parser.add_argument(
        "--insert_frequency",
        type=float,
        required=False,
        default=0.15,
        help="Minimum insert frequency for a site to be considered to have a variant. Default is 0.15"
    )
    parser.add_argument(
        "--delete_frequency",
        type=float,
        required=False,
        default=0.15,
        help="Minimum delete frequency for a site to be considered to have a variant. Default is 0.15"
    )
    parser.add_argument(
        "--min_coverage_threshold",
        type=int,
        required=False,
        default=5,
        help="Minimum delete frequency for a site to be considered to have a variant. Default is 5"
    )
    parser.add_argument(
        "--candidate_support_threshold",
        type=int,
        required=False,
        default=2,
        help="Minimum number of reads supporting a variant to be considered as a candidate. Default is 2"
    )
    parser.add_argument(
        "--candidate_frequency_threshold",
        type=float,
        required=False,
        default=0.10,
        help="Minimum frequency for a candidate to be considered to be a variant. Default is 0.10"
    )
    parser.add_argument(
        "--skip_indels",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs for inference. CUDA required. Default is False."
    )

    # Inference parameters
    parser.add_argument(
        "-bs",
        "--batch_size",
        type=int,
        required=False,
        default=128,
        help="Batch size for testing, default is 128. Suggested values: 256/512/1024. Default is 128"
    )
    parser.add_argument(
        "-g",
        "--gpu",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs for inference. CUDA required. Default is False."
    )
    parser.add_argument(
        "-per_gpu",
        "--callers_per_gpu",
        type=int,
        required=False,
        default=4,
        help="Number of callers to initialize per GPU. Default is 4"
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
        "--quantized",
        default=False,
        action='store_true',
        help="Use quantization for inference while on CPU inference mode. Speeds up inference. Default is False."
    )
    parser.add_argument(
        "-w",
        "--num_workers",
        type=int,
        required=False,
        default=8,
        help="Number of workers for loading images. Default is 8"
    )
    # find candidates arguments
    parser.add_argument(
        "--allowed_multiallelics",
        type=int,
        required=False,
        default=4,
        help="Number of maximum multialleleic variants allowed per site. Default is 4"
    )
    parser.add_argument(
        "--snp_p_value",
        required=False,
        type=float,
        default=0.1,
        help="Predicted value used for a SNP to be considered a candidate. Default is 0.4"
    )
    parser.add_argument(
        "--insert_p_value",
        required=False,
        type=float,
        default=0.3,
        help="Predicted value used for a insert to be considered a candidate. Default is 0.1"
    )
    parser.add_argument(
        "--delete_p_value",
        required=False,
        type=float,
        default=0.2,
        help="Predicted value used for a delete to be considered a candidate. Default is 0.2"
    )
    profile_group = parser.add_mutually_exclusive_group(required=True)
    profile_group.add_argument("--ont",
                               default=False,
                               action='store_true',
                               help="Set to call variants on Oxford Nanopore reads.")
    profile_group.add_argument("--hifi",
                               default=False,
                               action='store_true',
                               help="Set to call variants on PacBio HiFi reads.")
    profile_group.add_argument("--clr",
                               default=False,
                               action='store_true',
                               help="Set to call variants on PacBio CLR reads.")
    profile_group.add_argument("--custom",
                               default=False,
                               action='store_true',
                               help="Set to call variants with user-defined parameters.")
    return parser