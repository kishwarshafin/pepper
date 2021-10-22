def add_find_candidates_arguments(parser):
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
        help="Path to directory containing HDF files."
    )
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
        "-t",
        "--threads",
        required=True,
        type=int,
        help="Number of threads."
    )
    parser.add_argument(
        "-hp",
        "--use_hp_info",
        default=False,
        action='store_true',
        help="If set then haplotype-aware mode will be enabled."
    )
    parser.add_argument(
        "--allowed_multiallelics",
        type=int,
        required=False,
        default=None,
        help="Number of maximum multialleleic variants allowed per site."
    )
    parser.add_argument(
        "--snp_p_value",
        required=False,
        type=float,
        default=None,
        help="Predicted value used for a SNP to be considered a candidate."
    )
    parser.add_argument(
        "--insert_p_value",
        required=False,
        type=float,
        default=None,
        help="Predicted value used for a insert to be considered a candidate."
    )
    parser.add_argument(
        "--delete_p_value",
        required=False,
        type=float,
        default=None,
        help="Predicted value used for a delete to be considered a candidate."
    )
    parser.add_argument(
        "--freq_based",
        default=False,
        action='store_true',
        help="If set then frequency based variants will be reported."
    )
    parser.add_argument(
        "--freq",
        required=False,
        type=float,
        default=0.10,
        help="If frequency based variant finding in enabled then this frequency will be the threshold."
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
    return parser