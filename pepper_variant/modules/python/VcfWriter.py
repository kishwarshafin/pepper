import pysam
from pysam import VariantFile, VariantHeader
from pepper_variant.build import PEPPER_VARIANT
import collections
import math
import numpy as np
Candidate = collections.namedtuple('Candidate', 'chromosome_name pos_start pos_end ref '
                                                'alternate_alleles allele_depths '
                                                'allele_frequencies genotype qual gq predictions')


class VCFWriter:
    def __init__(self, all_contigs, reference_file_path, sample_name, output_dir, filename_full, filename_pepper, filename_variant_calling):
        self.fasta_handler = PEPPER_VARIANT.FASTA_handler(reference_file_path)
        contigs = self.fasta_handler.get_chromosome_names()
        # contigs = [x for x in contigs if x in all_contigs]
        self.contigs = contigs
        self.vcf_header = self.get_vcf_header(sample_name, contigs)
        self.output_dir = output_dir

        self.full_vcf_file_name = self.output_dir + filename_full + '.vcf.gz'
        self.pepper_vcf_file_name = self.output_dir + filename_pepper + '.vcf.gz'
        self.variant_vcf_file_name = self.output_dir + filename_variant_calling + '.vcf.gz'
        self.snp_variant_vcf_file_name = self.output_dir + filename_variant_calling + '_SNPs.vcf.gz'
        self.indel_variant_vcf_file_name = self.output_dir + filename_variant_calling + '_INDEL.vcf.gz'

        self.vcf_file_full = VariantFile(self.full_vcf_file_name, 'w', header=self.vcf_header)
        self.vcf_file_pepper = VariantFile(self.pepper_vcf_file_name, 'w', header=self.vcf_header)
        self.vcf_file_variant_calling = VariantFile(self.variant_vcf_file_name, 'w', header=self.vcf_header)
        self.vcf_file_variant_calling_snp = VariantFile(self.snp_variant_vcf_file_name, 'w', header=self.vcf_header)
        self.vcf_file_variant_calling_indel = VariantFile(self.indel_variant_vcf_file_name, 'w', header=self.vcf_header)

    def __del__(self):
        # close the files
        self.vcf_file_full.close()
        self.vcf_file_pepper.close()
        self.vcf_file_variant_calling.close()
        self.vcf_file_variant_calling_snp.close()
        self.vcf_file_variant_calling_indel.close()

        # index the files
        pysam.tabix_index(self.full_vcf_file_name, preset="vcf", force=True)
        pysam.tabix_index(self.pepper_vcf_file_name, preset="vcf", force=True)
        pysam.tabix_index(self.variant_vcf_file_name, preset="vcf", force=True)
        pysam.tabix_index(self.snp_variant_vcf_file_name, preset="vcf", force=True)
        pysam.tabix_index(self.indel_variant_vcf_file_name, preset="vcf", force=True)

    def candidate_list_to_variant(self, candidates, options):
        candidates = sorted(candidates, key=lambda x: (x[5], x[8]), reverse=True)
        if len(candidates) > options.allowed_multiallelics:
            candidates = candidates[:options.allowed_multiallelics]
        max_ref_length = 0
        max_ref_allele = ''

        # find the maximum reference allele
        for candidate in candidates:
            contig, ref_start, ref_end, ref_allele, alt_allele, genotype, depth, variant_allele_support, genotype_probability, predictions, non_alt_predictions, in_repeat = candidate
            if len(ref_allele) > max_ref_length:
                max_ref_length = len(ref_allele)
                max_ref_allele = ref_allele

        normalized_candidates = []
        for candidate in candidates:
            contig, ref_start, ref_end, ref_allele, alt_allele, genotype, depth, variant_allele_support, genotype_probability, predictions, non_alt_predictions, in_repeat = candidate
            suffix_needed = 0
            if len(ref_allele) < max_ref_length:
                suffix_needed = max_ref_length - len(ref_allele)

            if suffix_needed > 0:
                suffix_seq = max_ref_allele[-suffix_needed:]
                ref_allele = ref_allele + suffix_seq
                alt_allele = [alt + suffix_seq for alt in alt_allele]

            normalized_candidates.append((contig, ref_start, ref_end, ref_allele, alt_allele, genotype, depth, variant_allele_support, genotype_probability, predictions, non_alt_predictions, in_repeat))

        candidates = normalized_candidates
        gt_qual = -1.0
        genotype_hp1 = []
        genotype_hp2 = []

        all_initialized = False
        site_contig = ''
        site_ref_start = 0
        site_ref_end = 0
        site_ref_allele = ''
        site_depth = 0
        site_alts = []
        site_supports = []
        site_qualities = []
        site_in_repeat = False
        site_non_alt_predictions = []

        for i, candidate in enumerate(candidates):
            contig, ref_start, ref_end, ref_allele, alt_allele, genotype, depth, variant_allele_support, genotype_probability, predictions, non_alt_predictions, in_repeat = candidate
            # print(contig, ref_start, ref_end, ref_allele, alt_allele, genotype, depth, variant_allele_support, genotype_probability, predictions)
            site_in_repeat = in_repeat or site_in_repeat
            predicted_genotype = np.argmax(predictions)
            if predicted_genotype != 0:
                if gt_qual < 0:
                    gt_qual = predictions[predicted_genotype]
                else:
                    gt_qual = min(gt_qual, predictions[predicted_genotype])
            else:
                if gt_qual < 0:
                    # gt_qual = 1.0 - predictions[0]
                    gt_qual = max(predictions[1], predictions[2])

            if not all_initialized:
                site_contig = contig
                site_ref_start = ref_start
                site_ref_end = ref_start + len(ref_allele)
                site_ref_allele = ref_allele
                site_depth = depth
                all_initialized = True

            site_depth = min(site_depth, depth)
            site_alts.append(alt_allele[0])
            site_supports.append(variant_allele_support[0])
            site_qualities.append(genotype_probability)
            site_non_alt_predictions.extend(non_alt_predictions)

            # het
            if predicted_genotype == 1:
                genotype_hp1.append(i+1)
            # hom-alt
            elif predicted_genotype == 2:
                genotype_hp1.append(i+1)
                genotype_hp2.append(i+1)

        if 0 < len(genotype_hp1) + len(genotype_hp2) <= 2:
            gt = genotype_hp1 + genotype_hp2
            if len(gt) == 1:
                gt = [0, gt[0]]
        else:
            gt = [0, 0]

        # print(site_contig, site_ref_start, site_ref_end, site_ref_allele, site_alts, gt, site_depth, site_supports, gt_qual)
        return site_contig, site_ref_start, site_ref_end, site_ref_allele, site_alts, gt, site_depth, site_supports, gt_qual, site_non_alt_predictions, site_in_repeat

    def write_vcf_records(self, variants_list, options):
        total_variants, total_variants_pepper, total_variants_variant_calling, total_variants_variant_calling_snp, total_variants_variant_calling_indel = 0, 0, 0, 0, 0

        last_position = -1
        for contig, position in sorted(variants_list):
            all_candidates = variants_list[(contig, position)]

            contig, ref_start, ref_end, ref_seq, alleles, genotype, depth, variant_allele_support, genotype_probability, non_alt_predictions, site_in_repeat = self.candidate_list_to_variant(all_candidates, options)
            # print("RETURNED", contig, ref_start, ref_end, ref_seq, alleles, genotype, depth, variant_allele_support, genotype_probability)
            if len(alleles) <= 0:
                continue
            if ref_start == last_position:
                continue
            max_alt_len = max(len(ref_seq), max([len(x) for x in alleles]))
            last_position = ref_start
            alleles = tuple([ref_seq]) + tuple(alleles)
            qual = max(1, int(-10 * math.log10(max(0.000000001, 1.0 - genotype_probability))))
            alt_qual = max(1, int(-10 * math.log10(max(0.000000001, 1.0 - genotype_probability))))
            failed_variant = False
            is_snp = False
            if max_alt_len == 1:
                # this is SNP
                is_snp = True
                if not site_in_repeat and qual <= options.snp_q_cutoff:
                    failed_variant = True
                elif site_in_repeat and qual <= options.snp_q_cutoff_in_lc:
                    failed_variant = True
            else:
                if not site_in_repeat and qual <= options.indel_q_cutoff:
                    failed_variant = True
                elif site_in_repeat and qual <= options.indel_q_cutoff_in_lc:
                    failed_variant = True

            selected_for_variant_calling = False
            # Mode 1 are variants we will NOT re-genotype. Mode 2 is the variants selected for re-genotyping.
            if genotype == [0, 0] or failed_variant:
                # all variants that could not be genotyped, we will always re-genotype them
                selected_for_variant_calling = True

            vafs = [round(ad/max(1, depth), 3) for ad in variant_allele_support]

            # print(str(contig), ref_start, ref_end, qual, alleles, genotype, alt_qual, alt_qual, depth, variant_allele_support, vafs)
            if site_in_repeat:
                rep = "1"
            else:
                rep = "0"

            # always put things in all vcf
            if genotype == [0, 0]:
                vcf_record = self.vcf_file_full.new_record(contig=str(contig), start=ref_start,
                                                           stop=ref_end, id='.', qual=qual,
                                                           filter='refCall', alleles=alleles, GT=genotype,
                                                           AP=non_alt_predictions, GQ=alt_qual, DP=depth, AD=variant_allele_support, VAF=vafs,
                                                           REP=rep)
            else:
                vcf_record = self.vcf_file_full.new_record(contig=str(contig), start=ref_start,
                                                           stop=ref_end, id='.', qual=qual,
                                                           filter='PASS', alleles=alleles, GT=genotype,
                                                           AP=non_alt_predictions, GQ=qual, DP=depth, AD=variant_allele_support, VAF=vafs,
                                                           REP=rep)

            self.vcf_file_full.write(vcf_record)
            total_variants += 1

            if selected_for_variant_calling:
                if is_snp:
                    self.vcf_file_variant_calling_snp.write(vcf_record)
                    total_variants_variant_calling_snp += 1
                else:
                    self.vcf_file_variant_calling_indel.write(vcf_record)
                    total_variants_variant_calling_indel += 1

                self.vcf_file_variant_calling.write(vcf_record)
                total_variants_variant_calling += 1
            else:
                self.vcf_file_pepper.write(vcf_record)
                total_variants_pepper += 1

        return total_variants, total_variants_pepper, total_variants_variant_calling, total_variants_variant_calling_snp, total_variants_variant_calling_indel

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
                 ('Number', "A"),
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
        items = [('ID', "REP"),
                 ('Number', 1),
                 ('Type', 'String'),
                 ('Description', "If set to 1 then variant site is considered to be ina LowCompexity repeat region")]
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
