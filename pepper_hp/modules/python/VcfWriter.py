from pysam import VariantFile, VariantHeader
from pepper_hp.build import PEPPER_HP
from tqdm import tqdm
import collections
Candidate = collections.namedtuple('Candidate', 'chromosome_name pos_start pos_end ref '
                                                'alternate_alleles allele_depths '
                                                'allele_frequencies genotype qual gq predictions')


class VCFWriter:
    def __init__(self, reference_file_path, contigs, sample_name, output_dir, filename):
        self.fasta_handler = PEPPER_HP.FASTA_handler(reference_file_path)
        self.contigs = contigs
        self.vcf_header = self.get_vcf_header(sample_name, contigs)
        self.output_dir = output_dir
        self.filename = filename

    def write_vcf_records(self, variants_list):
        last_position = -1
        with VariantFile(self.output_dir + self.filename + '.vcf.gz', 'w', header=self.vcf_header) as vcf_file:
            for called_variant in tqdm(variants_list):
                contig, ref_start, ref_end, ref_seq, alleles, genotype, dps, gqs, ads, non_ref_prob = called_variant
                if ref_start == last_position:
                    continue
                last_position = ref_start
                alleles = tuple([ref_seq]) + tuple(alleles)
                # qual = -10 * math.log10(max(0.000001, 1.0 - max(0.0001, non_ref_prob)))
                qual = non_ref_prob

                # phred_gqs = []
                # for gq in gqs:
                #     phred_gq = -10 * math.log10(max(0.000001, 1.0 - max(0.0001, gq)))
                #     phred_gqs.append(phred_gq)
                vafs = [round(ad/max(1, max(dps)), 3) for ad in ads]
                if genotype == [0, 0]:
                    vcf_record = vcf_file.new_record(contig=str(contig), start=ref_start,
                                                     stop=ref_end, id='.', qual=qual,
                                                     filter='refCall', alleles=alleles, GT=genotype, GQ=min(gqs), VAF=vafs)
                else:
                    vcf_record = vcf_file.new_record(contig=str(contig), start=ref_start,
                                                     stop=ref_end, id='.', qual=qual,
                                                     filter='PASS', alleles=alleles, GT=genotype, GQ=min(gqs), VAF=vafs)
                vcf_file.write(vcf_record)

    def get_vcf_header(self, sample_name, contigs):
        header = VariantHeader()

        sqs = self.fasta_handler.get_chromosome_names()
        for sq in sqs:
            if sq not in contigs:
                continue
            sq_id = sq
            ln = self.fasta_handler.get_chromosome_sequence_length(sq)
            items = [('ID', sq_id),
                     ('length', ln)]
            items = [('ID', sq_id)]
            header.add_meta(key='contig', items=items)

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

        header.add_sample(sample_name)

        return header
