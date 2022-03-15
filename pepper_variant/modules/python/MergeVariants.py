import sys
from datetime import datetime
from collections import defaultdict
from pysam import VariantFile
from pepper_variant.modules.python.MergedVcfWriter import VCFWriter


def merge_vcf_records(options):
    pepper_vcf_file = VariantFile(options.vcf_pepper)
    all_pepper_records = pepper_vcf_file.fetch()
    positional_dv_records = defaultdict()
    total_records = 0
    output_vcf_name = 'PEPPER_MARGIN_DEEPVARIANT_OUTPUT.vcf.gz'

    if options.vcf_deepvariant:
        deepvariant_vcf_file = VariantFile(options.vcf_deepvariant)
        all_deepvariant_records = deepvariant_vcf_file.fetch()

        for record in all_deepvariant_records:
            positional_dv_records[(record.chrom, record.pos)] = record
            total_records += 1
        dv_samples = list(deepvariant_vcf_file.header.samples)
    else:
        deepvariant_vcf_file_snps = VariantFile(options.vcf_deepvariant_snps)
        deepvariant_vcf_file_indels = VariantFile(options.vcf_deepvariant_indels)
        all_snp_records = deepvariant_vcf_file_snps.fetch()
        all_indel_records = deepvariant_vcf_file_indels.fetch()

        # record all snps
        for record in all_snp_records:
            positional_dv_records[(record.chrom, record.pos)] = record
            total_records += 1
        # record all indels
        for record in all_indel_records:
            positional_dv_records[(record.chrom, record.pos)] = record
            total_records += 1
        dv_samples = list(deepvariant_vcf_file_snps.header.samples)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL VARIANTS IN DeepVariant: " + str(total_records) + "\n")

    # get all contigs from PEPPER
    contigs = []
    pepper_samples = list(pepper_vcf_file.header.samples)
    print(pepper_samples)
    print(dv_samples)

    if len(pepper_samples) > 1 or len(dv_samples) > 1:
        raise ValueError("ERROR: VARIANT FILE HAS MORE THAN ONE SAMPLE.")
    if pepper_samples[0] != dv_samples[0]:
        raise ValueError("ERROR: SAMPLE NAMES IN TWO VCFs DO NOT MATCH.")
    sample = pepper_samples[0]

    for x in pepper_vcf_file.header.records:
        if x.type == "CONTIG":
            contigs.append((x['ID'], x['length']))

    # output file
    vcf_out = VCFWriter(contigs, sample, options.output_dir, output_vcf_name)
    total_pepper_calls = 0
    total_dv_calls = 0
    total_pass_calls = 0
    for record in all_pepper_records:
        is_dv = False
        if (record.chrom, record.pos) in positional_dv_records:
            final_record = positional_dv_records[(record.chrom, record.pos)]
            is_dv = True
            total_dv_calls += 1
        else:
            final_record = record
            total_pepper_calls += 1

        if 'PASS' in final_record.filter.keys():
            total_pass_calls += 1

        vcf_out.write_vcf_records(final_record, sample, is_dv)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL VARIANTS FROM PEPPER: " + str(total_pepper_calls) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL VARIANTS FROM DEEPVARIANT: " + str(total_dv_calls) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL PASS VARIANTS: " + str(total_pass_calls) + "\n")
