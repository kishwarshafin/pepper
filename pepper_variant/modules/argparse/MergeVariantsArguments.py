def add_merge_variants_arguments(parser):
    """
    Add arguments to a parser for sub-command "stitch"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-vp",
        "--vcf_pepper",
        type=str,
        required=True,
        help="Path to VCF file from PEPPER SNP."
    )
    parser.add_argument(
        "-vd",
        "--vcf_deepvariant",
        type=str,
        required=True,
        help="Path to VCF file from DeepVariant."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory."
    )
    return parser