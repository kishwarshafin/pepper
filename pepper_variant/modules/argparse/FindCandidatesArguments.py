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
    return parser