from collections import defaultdict
from pysam import VariantFile


def merge_vcf_records(pepper_vcf_file, deepvariant_vcf_file, quality_threshold, output_directory):
    pepper_vcf_file = VariantFile(pepper_vcf_file)
    all_pepper_records = pepper_vcf_file.fetch()
    output_vcf = output_directory + "/" + 'PEPPER_DEEPVARIANT_MERGED.vcf.gz'

    positional_pepper_records = defaultdict(lambda: defaultdict())
    for record in all_pepper_records:
        if record.qual >= quality_threshold:
            positional_pepper_records[record.chrom][record.pos] = record

    deepvariant_vcf_file = VariantFile(deepvariant_vcf_file)
    all_deepvariant_records = deepvariant_vcf_file.fetch()
    vcf_out = VariantFile(output_vcf, 'w', header=deepvariant_vcf_file.header)

    genotype = (0, 0)
    for record in all_deepvariant_records:
        for sample_name, sample_items in record.samples.items():
            sample_items = sample_items.items()
            for name, value in sample_items:
                if name == 'GT':
                    genotype = value

        if 'PASS' not in record.filter.keys() or genotype == (0, 0):
            if record.pos in positional_pepper_records[record.chrom].keys():
                pepper_record = positional_pepper_records[record.chrom][record.pos]
                # Fetch PEPPER's genotype
                pepper_genotype = (0, 0)

                for sample in pepper_record.samples:
                    pepper_genotype = pepper_record.samples[sample]['GT']

                # apply PEPPER's genotype to DeepVariant and make it PASS
                if 'PASS' in pepper_record.filter.keys() and record.alleles == pepper_record.alleles:
                    for sample in pepper_record.samples:
                        record.samples[sample]['GT'] = pepper_genotype
                        record.qual = 5

                record.filter.add('PASS')

        vcf_out.write(record)
