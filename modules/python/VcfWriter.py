from pysam import VariantFile, VariantHeader
from modules.python.bam_handler import BamHandler
import time
import collections
Candidate = collections.namedtuple('Candidate', 'chromosome_name pos_start pos_end ref '
                                                'alternate_alleles allele_depths '
                                                'allele_frequencies genotype qual gq predictions')

class VCFWriter:
    def __init__(self, bam_file_path, sample_name, output_dir, contigs):
        self.bam_handler = BamHandler(bam_file_path)
        bam_file_name = bam_file_path.rstrip().split('/')[-1].split('.')[0]
        vcf_header = self.get_vcf_header(sample_name, contigs)
        time_str = time.strftime("%m%d%Y_%H%M%S")

        self.vcf_file = VariantFile(output_dir + bam_file_name + '_' + time_str + '.vcf', 'w', header=vcf_header)

    def add_variants(self, called_variants):
        record_set = set()
        called_variants = sorted(called_variants, key=lambda tup: tup[1])

        for variant in called_variants:
            chromosome_name, pos, pos_end, ref, alternate_alleles, genotype = variant

            if (chromosome_name, pos) in record_set:
                continue

            record_set.add((chromosome_name, pos))
            alleles = tuple([ref]) + tuple(alternate_alleles)
            qual = 20
            gq = 20
            vcf_record = self.vcf_file.new_record(contig=chromosome_name, start=pos,
                                                  stop=pos_end, id='.', qual=qual,
                                                  filter='PASS', alleles=alleles, GT=genotype, GQ=gq)
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
        bam_sqs = self.bam_handler.get_header_sq()
        for sq in bam_sqs:
            id = sq['SN']
            ln = sq['LN']
            if id not in contigs:
                continue
            items = [('ID', id),
                     ('length', ln)]
            header.add_meta(key='contig', items=items)

        header.add_sample(sample_name)

        return header
