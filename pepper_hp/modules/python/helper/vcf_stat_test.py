import argparse
from collections import defaultdict
from pysam import VariantFile
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def plot_distributions(vafs, true_vafs, false_vafs):
    sns.set(rc={"figure.figsize": (20, 10)})
    sns.set_style("white")
    # plt.hist(vafs, bins=100, color='blue', alpha=0.4)
    plt.hist([true_vafs, false_vafs], bins=100, density=False, histtype='bar', color=['green', 'red'], alpha=0.4, stacked=True, label=['True variants', 'False positives'])
    # plt.hist(false_vafs, bins=100, color='red', alpha=0.4, stacked=True)
    plt.xlim((0.00, 1.10))
    plt.ylim((0, 5000))
    plt.legend(fontsize='x-large')
    # plt.xticks(np.arange(0, 1, step=0.10), fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Allele frequency", fontsize='24')
    plt.ylabel("Count", fontsize='24')

    plt.title("Stacked histogram showing TP and FP distribution at different frequency intervals.", fontsize='20')
    # plt.show()
    # exit()
    output_file_name = "./VAF_distribution.png"
    plt.savefig(output_file_name, format='png', dpi=300, quality=95)

def plot_distributions_allele_probs(all_aps, true_aps, false_aps):
    sns.set(rc={"figure.figsize": (20, 10)})
    sns.set_style("white")
    # plt.hist(vafs, bins=100, color='blue', alpha=0.4)
    # plt.hist([true_aps, false_aps], bins=100, density=False, histtype='bar', color=['green', 'red'], alpha=0.4, stacked=True, label=['True variants', 'False positives'])
    plt.hist(false_aps, color='red', bins=1000, alpha=0.4)
    plt.hist(true_aps, color='green', bins=1000, alpha=0.5)
    # plt.xlim((0.00, 0.20))
    # plt.ylim((0, 200))
    plt.legend(fontsize='x-large')
    # plt.xticks(np.arange(0, 1, step=0.10), fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Allele probability", fontsize='24')
    plt.ylabel("Count", fontsize='24')

    plt.title("Stacked histogram showing TP and FP distribution at different frequency intervals.", fontsize='20')
    plt.show()
    exit()
    output_file_name = "./VAF_distribution.png"
    plt.savefig(output_file_name, format='png', dpi=300, quality=95)


def plot_distributions_non_ref(true_qvals, false_qvals):
    sns.set(rc={"figure.figsize": (20, 10)})
    sns.set_style("white")
    # plt.hist(vafs, bins=100, color='blue', alpha=0.4)
    plt.hist([true_qvals, false_qvals], bins=100, density=False, histtype='bar', color=['green', 'red'], alpha=0.4, stacked=True, label=['True variants', 'False positives'])
    # plt.hist(false_vafs, bins=100, color='red', alpha=0.4, stacked=True)
    # plt.xlim((0.00, 0.2))
    # plt.ylim((0, 500))
    plt.legend(fontsize='x-large')
    # plt.xticks(np.arange(0, 1, step=0.10), fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("P-value", fontsize='24')
    plt.ylabel("Count", fontsize='24')

    plt.title("Stacked histogram showing TP and FP distribution at different P-value intervals.", fontsize='20')
    plt.show()
    exit()
    output_file_name = "./P_val_distribution.png"
    plt.savefig(output_file_name, format='png', dpi=300, quality=95)


def plot_vaf_and_q_2d(true_allele_q_vaf, false_allele_q_vaf, q_vaf):
    sns.set(rc={"figure.figsize": (20, 10)})
    sns.set_style("white")
    q, vaf, pred = zip(*q_vaf)
    plt.scatter(vaf, q, c=pred)
    plt.xlabel("VAF", fontsize='24')
    plt.ylabel("Q-value", fontsize='24')
    plt.xlim((0.00, 0.25))
    plt.ylim((0.00, 0.025))
    plt.show()


def calculate_stats(truth_vcf, vcf):
    untagged_vcf = VariantFile(vcf)
    positional_norm_q = defaultdict(lambda: defaultdict(float))
    positional_aps = defaultdict(lambda: defaultdict(list))
    for rec in untagged_vcf.fetch("chr20"):
        positional_norm_q[rec.chrom][rec.pos] = rec.qual
        alts = rec.alts

        # get the allele probability from the original vcf
        for sample in rec.samples:
            if 'AP' in rec.samples[sample].keys():
                aps = list(rec.samples[sample]['AP'])
            else:
                aps = [0.0] * len(alts)
            positional_aps[rec.chrom][rec.pos] = aps

    vcf_in1 = VariantFile(truth_vcf)

    vafs_of_true_alleles = list()
    vafs_of_false_alleles = list()
    non_ref_of_true_alleles = list()
    non_ref_of_false_alleles = list()

    ap_of_true_alleles = list()
    ap_of_false_alleles = list()
    total_alts = 0
    total_true_calls = 0
    total_false_calls = 0
    true_allele_q_vaf = list()
    false_allele_q_vaf = list()

    ap_vaf_list = list()

    VAF_Threshold = 1.10

    all_allele_frequencies = list()
    all_aps = list()

    total_recs = 0
    total_multi_allelic_sites = 0
    print("CONTIG" + "\t" + "POS" + "\t" + "REF" + "\t" + "ALT" + "\t" + "ALLELE_FREQ" + "\t" + "ALLELE_WEIGHT" + "\t" + "NON_REF_PROB" + "\t" + "T/F")
    for rec in vcf_in1.fetch("chr20"):
        total_recs += 1
        # ['__class__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__',
        # 'alleles', 'alts', 'chrom', 'contig', 'copy', 'filter', 'format', 'header', 'id', 'info', 'pos', 'qual', 'ref', 'rid', 'rlen', 'samples', 'start', 'stop', 'translate']
        alts = rec.alts
        is_indel = False
        if len(rec.ref) > 1:
            is_indel = True
        for alt in alts:
            if len(alt) > 1:
                is_indel = True

        # toggle this to switch between switch and indels
        # if not is_indel:
        #     continue

        total_alts += len(alts)
        for sample in rec.samples:
            # allele frequencies
            if 'VAF' in rec.samples[sample].keys():
                vafs = list(rec.samples[sample]['VAF'])
            else:
                vafs = [0.0] * len(alts)

            # allele weights
            aps = positional_aps[rec.chrom][rec.pos]
            gts = list(rec.samples[sample]['GT'])

            true_index = []
            for gt in gts:
                if gt != 0:
                    true_index.append(gt-1)

            # generate the non-ref probability distribution
            if true_index:
                non_ref_of_true_alleles.append(positional_norm_q[rec.chrom][rec.pos])
            else:
                non_ref_of_false_alleles.append(positional_norm_q[rec.chrom][rec.pos])

            if len(alts) > 1:
                total_multi_allelic_sites += 1

            # generate allele probability distribution

            for i, (alt, vaf, ap) in enumerate(zip(alts, vafs, aps)):
                if i in true_index:
                    # this is a true allele tagged by the truth
                    print(rec.chrom + "\t" + str(rec.pos) + "\t" + rec.ref + "\t" + alt + "\t" + str(vaf) + "\t" + str(ap) + "\t" + str(positional_norm_q[rec.chrom][rec.pos]) + "\t" + "1")
                    vafs_of_true_alleles.append(vaf)
                    true_allele_q_vaf.append((ap, min(1.1, vaf)))

                    if vaf <= VAF_Threshold:
                        ap_vaf_list.append((ap, min(1.1, vaf), 'Green'))

                    ap_of_true_alleles.append(ap)

                    total_true_calls += 1
                else:
                    print(rec.chrom + "\t" + str(rec.pos) + "\t" + rec.ref + "\t" + alt + "\t" + str(vaf) + "\t" + str(ap) + "\t" + str(positional_norm_q[rec.chrom][rec.pos]) + "\t" + "0")
                    vafs_of_false_alleles.append(vaf)
                    false_allele_q_vaf.append((ap, min(1.1, vaf)))

                    if vaf <= VAF_Threshold:
                        ap_vaf_list.append((ap, min(1.1, vaf), 'Red'))

                    ap_of_false_alleles.append(ap)
                    total_false_calls += 1

            for vaf in vafs:
                all_allele_frequencies.append(round(vaf, 3))

            for ap in aps:
                all_aps.append(round(ap, 3))

    print("Records:\t", total_recs)
    print("Multi-alleleic:\t", total_false_calls, "(" + str(int(100 * (total_false_calls/total_alts))) + "%)")
    print("Alt alleles:\t", total_alts)
    print("True alleles:\t", total_true_calls, "(" + str(int(100 * (total_true_calls/total_alts))) + "%)")
    print("False alleles:\t", total_false_calls, "(" + str(int(100 * (total_false_calls/total_alts))) + "%)")

    # plot_distributions(all_allele_frequencies, vafs_of_true_alleles, vafs_of_false_alleles)
    # plot_distributions_allele_probs(all_aps, ap_of_true_alleles, ap_of_false_alleles)
    # plot_distributions_non_ref(non_ref_of_true_alleles, non_ref_of_false_alleles)
    # plot_vaf_and_q_2d(true_allele_q_vaf, false_allele_q_vaf, ap_vaf_list)



def add_merge_vcf_arguments(parser):
    parser.add_argument(
        "-vt",
        "--truth_vcf",
        type=str,
        required=True,
        help="VCF of haplotype 1."
    )
    parser.add_argument(
        "-v",
        "--vcf",
        type=str,
        required=True,
        help="VCF of haplotype 1."
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
    calculate_stats(FLAGS.truth_vcf, FLAGS.vcf)