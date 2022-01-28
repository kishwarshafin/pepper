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
        default=None,
        help="Minimum mapping quality for read to be considered valid. Default is 5"
    )
    parser.add_argument(
        "--min_snp_baseq",
        type=int,
        required=False,
        default=None,
        help="Minimum base quality for base to be considered valid for snp."
    )
    parser.add_argument(
        "--min_indel_baseq",
        type=int,
        required=False,
        default=None,
        help="Minimum base quality for base to be considered valid for indels."
    )
    parser.add_argument(
        "--snp_frequency",
        type=float,
        required=False,
        default=None,
        help="Minimum SNP frequency for a site to be considered to have a variant. Default is 0.10"
    )
    parser.add_argument(
        "--insert_frequency",
        type=float,
        required=False,
        default=None,
        help="Minimum insert frequency for a site to be considered to have a variant. Default is 0.15"
    )
    parser.add_argument(
        "--delete_frequency",
        type=float,
        required=False,
        default=None,
        help="Minimum delete frequency for a site to be considered to have a variant. Default is 0.15"
    )
    parser.add_argument(
        "--min_coverage_threshold",
        type=int,
        required=False,
        default=None,
        help="Minimum delete frequency for a site to be considered to have a variant. Default is 5"
    )
    parser.add_argument(
        "--candidate_support_threshold",
        type=int,
        required=False,
        default=None,
        help="Minimum number of reads supporting a variant to be considered as a candidate. Default is 2"
    )
    parser.add_argument(
        "--snp_candidate_frequency_threshold",
        type=float,
        required=False,
        default=None,
        help="Minimum frequency for a SNP candidate to be considered to be a variant."
    )
    parser.add_argument(
        "--indel_candidate_frequency_threshold",
        type=float,
        required=False,
        default=None,
        help="Minimum frequency for a SNP candidate to be considered to be a variant."
    )
    parser.add_argument(
        "--skip_indels",
        default=False,
        action='store_true',
        help="If set then INDEL calling will be skipped."
    )
    profile_group = parser.add_mutually_exclusive_group(required=True)
    profile_group.add_argument("--ont_r9_guppy5_sup",
                               default=False,
                               action='store_true',
                               help="Set to call variants on R9.4.1 Guppy SUP caller Oxford Nanopore reads.")
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