import argparse
from collections import defaultdict
from pysam import VariantFile


def merg_vcf(h1_vcf, h2_vcf, output_dir, merge_genotype):

    vcf_positional_dict = defaultdict(lambda: defaultdict(list))
    vcf_in1 = VariantFile(h1_vcf)
    vcf_out = VariantFile(output_dir + 'merged_file.vcf', 'w', header=vcf_in1.header)

    for rec in vcf_in1.fetch():
        # ['__class__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__',
        # 'alleles', 'alts', 'chrom', 'contig', 'copy', 'filter', 'format', 'header', 'id', 'info', 'pos', 'qual', 'ref', 'rid', 'rlen', 'samples', 'start', 'stop', 'translate']
        if 'PASS' in rec.filter.keys():
            vcf_positional_dict[rec.chrom][rec.pos].append(rec)

    vcf_in2 = VariantFile(h2_vcf)
    for rec in vcf_in2.fetch():
        if 'PASS' in rec.filter.keys():
            vcf_positional_dict[rec.chrom][rec.pos].append(rec)

    for chrom in vcf_positional_dict.keys():
        for pos in sorted(vcf_positional_dict[chrom].keys()):
            # this means that merging is needed at this position
            if len(vcf_positional_dict[chrom][pos]) == 1:
                for var in vcf_positional_dict[chrom][pos]:
                    vcf_out.write(var)
            elif len(vcf_positional_dict[chrom][pos]) > 1:
                longest_ref = vcf_positional_dict[chrom][pos][0].ref
                longest_var = vcf_positional_dict[chrom][pos][0]
                for var in vcf_positional_dict[chrom][pos]:
                    if len(var.ref) > len(longest_ref):
                        longest_ref = var.ref
                        longest_var = var

                alts = [longest_ref]
                gq = -1.0
                qual = -1.0
                gts = []
                for var in vcf_positional_dict[chrom][pos]:
                    for sample in var.samples:
                        if gq < 0:
                            gq = var.samples[sample]['GQ']
                        gq = min(gq, var.samples[sample]['GQ'])
                        if var.samples[sample]['GT'] != [0, 0]:
                            gts.append(var.samples[sample]['GT'])

                    var_alts = list(var.alts)
                    var_ref = var.ref
                    if qual < 0:
                        qual = var.qual
                    qual = min(qual, var.qual)

                    ref_suffix = longest_ref[len(var_ref):]
                    for alt in var_alts:
                        if alt + ref_suffix not in alts and len(alt+ref_suffix) > 0:
                            alts.append(alt + ref_suffix)

                if len(alts) == 2:
                    if merge_genotype:
                        if len(gts) == 2:
                            genotype = [1, 1]
                        else:
                            genotype = [0, 1]
                    else:
                        genotype = gts[0]
                else:
                    genotype = [1, 2]

                vcf_record = vcf_out.new_record(contig=longest_var.contig,
                                                start=longest_var.start,
                                                stop=longest_var.stop,
                                                id=longest_var.id,
                                                qual=qual,
                                                filter=longest_var.filter,
                                                alleles=alts,
                                                GT=genotype,
                                                GQ=gq)
                vcf_out.write(vcf_record)


def add_merge_vcf_arguments(parser):
    parser.add_argument(
        "-v1",
        "--vcf_h1",
        type=str,
        required=True,
        help="VCF of haplotype 1."
    )
    parser.add_argument(
        "-v2",
        "--vcf_h2",
        type=str,
        required=True,
        help="VCF of haplotype 1."
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        required=True,
        help="FASTA file containing the reference assembly."
    )
    parser.add_argument(
        "-m",
        "--merge_genotype",
        default=False,
        action='store_true',
        help="If true then this will treat two VCFs from two haplotypes."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory, if it doesn't exist it will be created."
    )
    return parser


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="PEPPER is a RNN based polisher for polishing ONT-based assemblies. "
                                                 "It works in three steps:\n"
                                                 "1) make_images: This module takes alignment file and coverts them"
                                                 "to HDF5 files containing summary statistics.\n"
                                                 "2) run_inference: This module takes the summary images and a"
                                                 "trained neural network and generates predictions per base.\n"
                                                 "3) find_snps: This module takes the inference files as input and "
                                                 "finds possible SNP sites.\n",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )
    add_merge_vcf_arguments(parser)
    FLAGS, unparsed = parser.parse_known_args()
    merg_vcf(FLAGS.vcf_h1, FLAGS.vcf_h2, FLAGS.output_dir, FLAGS.merge_genotype)