from pysam import VariantFile, VariantHeader
from pepper_variant.build import PEPPER_VARIANT
import collections
import math
Candidate = collections.namedtuple('Candidate', 'chromosome_name pos_start pos_end ref '
                                                'alternate_alleles allele_depths '
                                                'allele_frequencies genotype qual gq predictions')


class VCFWriter:
    def __init__(self, reference_file_path, contigs, sample_name, output_dir, filename):
        self.fasta_handler = PEPPER_VARIANT.FASTA_handler(reference_file_path)
        self.contigs = contigs
        self.vcf_header = self.get_vcf_header(sample_name, contigs)
        self.output_dir = output_dir
        self.filename = filename
        self.vcf_file = VariantFile(self.output_dir + self.filename + '.vcf', 'w', header=self.vcf_header)

    def write_vcf_records(self, variants_list):
        last_position = -1
        for called_variant in variants_list:
            contig, ref_start, ref_end, ref_seq, alleles, genotype, dps, alt_probs, alt_prob_h1s, alt_prob_h2s, non_ref_probs, ads, overall_non_ref_prob = called_variant

            if len(alleles) <= 0:
                continue
            if ref_start == last_position:
                continue
            last_position = ref_start
            alleles = tuple([ref_seq]) + tuple(alleles)
            qual = max(1, int(-10 * math.log10(max(0.000001, 1.0 - max(0.0001, overall_non_ref_prob)))))
            overall_non_ref_prob = qual

            vafs = [round(ad/max(1, max(dps)), 3) for ad in ads]
            if genotype == [0, 0]:
                vcf_record = self.vcf_file.new_record(contig=str(contig), start=ref_start,
                                                      stop=ref_end, id='.', qual=qual,
                                                      filter='refCall', alleles=alleles, GT=genotype,
                                                      AP=alt_probs, APM=max(alt_probs), AP1=alt_prob_h1s, AP2=alt_prob_h2s,
                                                      NR=non_ref_probs, NRM=max(non_ref_probs),
                                                      GQ=overall_non_ref_prob, VAF=vafs)
            else:
                vcf_record = self.vcf_file.new_record(contig=str(contig), start=ref_start,
                                                      stop=ref_end, id='.', qual=qual,
                                                      filter='PASS', alleles=alleles, GT=genotype,
                                                      AP=alt_probs, APM=max(alt_probs), AP1=alt_prob_h1s, AP2=alt_prob_h2s,
                                                      NR=non_ref_probs, NRM=max(non_ref_probs),
                                                      GQ=overall_non_ref_prob, VAF=vafs)

            self.vcf_file.write(vcf_record)

    def get_vcf_header(self, sample_name, contigs):
        header = VariantHeader()

        items = [('ID', "PASS"),
                 ('Description', "All filters passed")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "refCall"),
                 ('Description', "Call is homozygous")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "lowGQ"),
                 ('Description', "Low genotype quality")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "lowQUAL"),
                 ('Description', "Low variant call quality")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "conflictPos"),
                 ('Description', "Overlapping record")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "GT"),
                 ('Number', 1),
                 ('Type', 'String'),
                 ('Description', "Genotype")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "DP"),
                 ('Number', 1),
                 ('Type', 'Integer'),
                 ('Description', "Depth")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "AD"),
                 ('Number', 1),
                 ('Type', 'String'),
                 ('Description', "Depth")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "VAF"),
                 ('Number', "A"),
                 ('Type', 'Float'),
                 ('Description', "Variant allele fractions.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "AP"),
                 ('Number', "A"),
                 ('Type', 'Float'),
                 ('Description', "Maximum variant allele probability for each allele.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "APM"),
                 ('Number', 1),
                 ('Type', 'Float'),
                 ('Description', "Maximum variant allele probability for the site.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "AP1"),
                 ('Number', "A"),
                 ('Type', 'Float'),
                 ('Description', "Maximum variant allele probability from hp1.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "AP2"),
                 ('Number', "A"),
                 ('Type', 'Float'),
                 ('Description', "Maximum variant allele probability from hp2.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "NR"),
                 ('Number', "A"),
                 ('Type', 'Float'),
                 ('Description', "Max probability of observing a non-ref allele for each allele.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "NRM"),
                 ('Number', 1),
                 ('Type', 'Float'),
                 ('Description', "Max probability of observing a non-ref allele for the site.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "GT"),
                 ('Number', 1),
                 ('Type', 'String'),
                 ('Description', "Genotype")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "GQ"),
                 ('Number', 1),
                 ('Type', 'Float'),
                 ('Description', "Genotype Quality")]
        header.add_meta(key='FORMAT', items=items)

        sqs = self.fasta_handler.get_chromosome_names()
        for sq in sqs:
            if sq not in contigs:
                continue
            sq_id = sq
            ln = self.fasta_handler.get_chromosome_sequence_length(sq)
            header.contigs.add(sq_id, length=ln)

        header.add_sample(sample_name)


        return header
