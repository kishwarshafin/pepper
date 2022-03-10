from pysam import VariantFile, VariantHeader
import pysam


class VCFWriter:
    def __init__(self, all_contigs, sample_name, output_dir, filename):
        self.vcf_header = self.get_vcf_header(sample_name, all_contigs)
        self.output_dir = output_dir
        self.filename = self.output_dir + "/" + filename
        self.vcf_file = VariantFile(self.filename, 'w', header=self.vcf_header)

    def __del__(self):
        # close the files
        self.vcf_file.close()

        # index the files
        pysam.tabix_index(self.filename, preset="vcf", force=True)

    def write_vcf_records(self, vcf_record, sample, is_deepvariant_call):
        if 'PASS' not in vcf_record.filter.keys():
            record_filter = 'refCall'
        else:
            record_filter = 'PASS'

        depth = vcf_record.samples[sample]['DP']
        genotype_quality = vcf_record.samples[sample]['GQ']
        genotype = vcf_record.samples[sample]['GT']
        variant_allele_frequencies = vcf_record.samples[sample]['VAF']
        if is_deepvariant_call:
            caller = 'DV'
            allele_depths = vcf_record.samples[sample]['AD'][1:]

        else:
            caller = 'P'
            allele_depths = vcf_record.samples[sample]['AD']

        vcf_record = self.vcf_file.new_record(contig=vcf_record.contig,
                                              start=vcf_record.start,
                                              stop=vcf_record.stop,
                                              id=vcf_record.id,
                                              qual=vcf_record.qual,
                                              filter=record_filter,
                                              alleles=vcf_record.alleles,
                                              GT=genotype,
                                              GQ=genotype_quality,
                                              DP=depth,
                                              AD=allele_depths,
                                              VAF=variant_allele_frequencies,
                                              C=caller)

        self.vcf_file.write(vcf_record)

    @staticmethod
    def get_vcf_header(sample_name, contigs):
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
                 ('Number', "A"),
                 ('Type', 'Integer'),
                 ('Description', "Allele depth")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "VAF"),
                 ('Number', "A"),
                 ('Type', 'Float'),
                 ('Description', "Variant allele fractions.")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "AP"),
                 ('Number', 1),
                 ('Type', 'Float'),
                 ('Description', "Maximum variant allele probability for each allele.")]
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
        items = [('ID', "C"),
                 ('Number', 1),
                 ('Type', 'String'),
                 ('Description', "Genotyper P=PEPPER DV=DeepVariant")]
        header.add_meta(key='FORMAT', items=items)

        for contig, length in contigs:
            header.contigs.add(contig, length=length)

        header.add_sample(sample_name)

        return header
