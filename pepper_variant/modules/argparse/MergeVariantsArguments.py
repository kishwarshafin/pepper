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
        required=False,
        help="Path to VCF file from DeepVariant."
    )
    parser.add_argument(
        "-vds",
        "--vcf_deepvariant_snps",
        type=str,
        required=False,
        help="Path to VCF file from DeepVariant with SNP records."
    )
    parser.add_argument(
        "-vdi",
        "--vcf_deepvariant_indels",
        type=str,
        required=False,
        help="Path to VCF file from DeepVariant with INDEL records."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory."
    )
    return parser