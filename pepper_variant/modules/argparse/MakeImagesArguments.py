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
        help="Number of threads to use."
    )
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

    return parser