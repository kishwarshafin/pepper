import argparse
from collections import defaultdict
from pysam import VariantFile
import re

def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def get_phase_block_length(vcf_file):

    vcf_PS_start_pos_dict = defaultdict(lambda: defaultdict(list))
    vcf_PS_end_pos_dict = defaultdict(lambda: defaultdict(list))
    chromosome_ps_blocks = defaultdict(list)
    vcf_in = VariantFile(vcf_file)

    for var in vcf_in.fetch():
        # ['__class__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__',
        # 'alleles', 'alts', 'chrom', 'contig', 'copy', 'filter', 'format', 'header', 'id', 'info', 'pos', 'qual', 'ref', 'rid', 'rlen', 'samples', 'start', 'stop', 'translate']
        if 'PASS' in var.filter.keys():
            for sample in var.samples:
                if 'PS' in var.samples[sample].keys():
                    if var.samples[sample]['PS'] not in vcf_PS_start_pos_dict[var.chrom].keys():
                        vcf_PS_start_pos_dict[var.chrom][var.samples[sample]['PS']] = var.start
                        vcf_PS_end_pos_dict[var.chrom][var.samples[sample]['PS']] = var.stop
                        chromosome_ps_blocks[var.chrom].append(var.samples[sample]['PS'])

                    vcf_PS_start_pos_dict[var.chrom][var.samples[sample]['PS']] = min(var.start, vcf_PS_start_pos_dict[var.chrom][var.samples[sample]['PS']])
                    vcf_PS_end_pos_dict[var.chrom][var.samples[sample]['PS']] = max(var.stop, vcf_PS_end_pos_dict[var.chrom][var.samples[sample]['PS']])

    contigs = sorted(chromosome_ps_blocks.keys(), key=natural_key)
    for contig in contigs:
        for ps_block in chromosome_ps_blocks[contig]:
            print(str(contig) + "," + str(ps_block) + "," + str(vcf_PS_start_pos_dict[contig][ps_block]) + "," + str(vcf_PS_end_pos_dict[contig][ps_block]))


def add_merge_vcf_arguments(parser):
    parser.add_argument(
        "-v1",
        "--vcf",
        type=str,
        required=True,
        help="VCF of haplotype 1."
    )
    return parser


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="HELPER SCRIPT\n",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )
    add_merge_vcf_arguments(parser)
    FLAGS, unparsed = parser.parse_known_args()
    get_phase_block_length(FLAGS.vcf)
