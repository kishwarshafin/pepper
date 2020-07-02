from pysam import VariantFile, VariantHeader
import time
import collections
from pepper_snp.build import PEPPER_SNP
Candidate = collections.namedtuple('Candidate', 'chromosome_name pos_start pos_end ref '
                                                'alternate_alleles allele_depths '
                                                'allele_frequencies genotype qual gq predictions')


class VCFWriter:
    def __init__(self, reference_path, sample_name, output_dir, contigs):
        self.fasta_handler = PEPPER_SNP.FASTA_handler(reference_path)
        vcf_header = self.get_vcf_header(sample_name, contigs)
        time_str = time.strftime("%m%d%Y_%H%M%S")

        self.vcf_file = VariantFile(output_dir + "CANDIDATES_PEPPER" + '_' + time_str + '.vcf', 'w', header=vcf_header)

    def get_genotype(self, ref, alt1, alt2):
        alt1_gt = 1
        alt2_gt = 2
        if ref == alt1 or alt1 == '*':
            alt1_gt = 0
        if ref == alt2 or alt2 == '*':
            alt2_gt = 0

        if alt1 == alt2:
            alt2_gt = alt1_gt
        gt = sorted([alt1_gt, alt2_gt])

        if gt == [0, 0]:
            return ref, [], [0, 0]
        if gt == [0, 1]:
            return ref, [alt1], [0, 1]
        if gt == [1, 1]:
            return ref, [alt1], [1, 1]
        if gt == [0, 2]:
            return ref, [alt2], [0, 1]
        if gt == [2, 2]:
            return ref, [alt2], [1, 1]
        if gt == [1, 2]:
            return ref, [alt1, alt2], [1, 2]

        return sorted([alt1_gt, alt2_gt])

    def get_alleles(self, ref_base, alt_predictions):
        alts1 = set()
        alts2 = set()
        for alt1, alt2 in alt_predictions:
            if alt1 != '*' and alt1 != ref_base:
                alts1.add(alt1)
            if alt2 != '*' and alt2 != ref_base:
                alts2.add(alt2)

        return list(alts1), list(alts2)

    def write_vcf_records(self, chromosome_name, called_variants, reference_dict, positions):
        for pos in sorted(positions):
            ref_base = reference_dict[pos]

            if ref_base == 'n' or ref_base == 'N':
                continue

            alts1, alts2 = self.get_alleles(ref_base, called_variants[pos])
            if alts1:
                alt1 = alts1[0]
            else:
                alt1 = ref_base

            if alts2:
                alt2 = alts2[0]
            else:
                alt2 = ref_base

            ref, alt_alleles, gt = self.get_genotype(ref_base, alt1, alt2)

            if gt == [0, 0]:
                continue
            # add extra alleles not used here
            for i in range(1, len(alts1)):
                alt_alleles.append(alts1[i])
            # add extra alleles not used
            for i in range(1, len(alts2)):
                alt_alleles.append(alts2[i])

            alleles = tuple([ref]) + tuple(set(alt_alleles))
            # print(str(chrm), st_pos, end_pos, qual, rec_filter, alleles, genotype, gq)
            vcf_record = self.vcf_file.new_record(contig=str(chromosome_name), start=pos,
                                                  stop=pos+1, id='.', qual=60,
                                                  filter='PASS', alleles=alleles, GT=gt, GQ=60)
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
