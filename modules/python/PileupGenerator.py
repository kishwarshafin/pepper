from build import HELEN
import numpy as np
from modules.python.Options import CandidateFinderOptions

SNP_CANDIDATE, IN_CANDIDATE, DEL_CANDIDATE = 1, 2, 3
# Genotype codes
HOM, HET, HOM_ALT = 0, 1, 2


class PileupGenerator:
    def __init__(self, fasta_handler, contig, start, end):
        self.fasta_handler = fasta_handler
        self.chromosome_name = contig
        self.region_start = start
        self.region_end = end

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    def get_label_of_allele(self, candidate, positional_vcf):
        """
        Given positional VCFs (IN, DEL, SNP), variant type and a candidate allele, return the try genotype.
        :param positional_vcf: Three dictionaries for each position
        :param candidate_allele: Candidate allele
        :param allele_types: Alt allele type: IN,DEL,SNP
        :return: genotype
        """
        gts = [0, 0]
        if candidate.pos not in positional_vcf:
            return gts

        for rec in positional_vcf[candidate.pos]:
            for alt in rec.alt_allele:
                # alt_ref = alt.ref
                # alt_alt = alt.alt_allele
                # if alt.alt_type == IN_CANDIDATE:
                #     alt_ref, alt_alt = self._resolve_suffix_for_insert(alt.ref, alt.alt_allele)
                # if alt.alt_type == DEL_CANDIDATE:
                #     alt_ref, alt_alt = self._resolve_suffix_for_delete(alt.ref, alt.alt_allele)

                # print(alt_ref, candidate.ref)
                # print(alt_alt, candidate.alt1)
                # print(alt.alt_type, candidate.alt1_type)
                # print(".....")
                if alt.ref == candidate.ref and \
                        alt.alt_allele == candidate.alt1 and \
                        alt.alt_type == candidate.alt1_type:
                    if rec.genotype[0] == 1 and rec.genotype[1] == 1:
                        gts[0] = HOM_ALT
                    elif rec.genotype[0] == 1 or rec.genotype[1] == 1:
                        gts[0] = HET
                # print(alt_ref, candidate.ref)
                # print(alt_alt, candidate.alt2)
                # print(alt.alt_type, candidate.alt2_type)
                if alt.ref == candidate.ref and \
                        alt.alt_allele == candidate.alt2 and \
                        alt.alt_type == candidate.alt2_type:
                    if rec.genotype[0] == 2 and rec.genotype[1] == 2:
                        gts[1] = HOM_ALT
                    elif rec.genotype[0] == 2 or rec.genotype[1] == 2:
                        gts[1] = HET
                # print("-------")
        return gts

    @staticmethod
    def get_combined_gt(gt):
        """
        Given two genotypes get the combined genotype. This is used to create labels for the third image.
        If two alleles have two different genotypes then the third genotype is inferred using this method.

        - If genotype1 is HOM then genotype of third image is genotype2
        - If genotype2 is HOM then genotype of third image is genotype1
        - If both gt are  HOM then genotype of third image is HOM
        - If genotype1, genotype2 both are HET then genotype of third image is HOM_ALT
        - If none of these cases match then we have an invalid genotype
        :param gt1: Genotype of first allele
        :param gt2: Genotype of second allele
        :return: genotype of image where both alleles are used together
        """
        gt1, gt2 = gt
        if gt1 == HOM and gt2 == HOM:
            return 0 # 0/0
        if gt1 == HET and gt2 == HOM:
            return 1 # 1/0
        if gt1 == HOM_ALT and gt2 == HOM:
            return 2 # 1/1
        if gt1 == HOM and gt2 == HET:
            return 3 # 0/2
        if gt2 == HOM_ALT and gt2 == HOM:
            return 4 # 2/2
        if gt1 == HET and gt2 == HET:
            return 5 # 1/2
        else:
            import sys
            sys.stderr.write("ERROR: GENOTYPES DIDNOT MATCH: ", gt1, gt2)
        return None

    def decode_image_row(self, read_array):
        global_base_color_reverse = {250: 'A', 30: 'C', 180: 'G', 100: 'T', 0: ' ', 10: 'N'}
        for pixel in read_array:
            base = pixel[0]
            print(global_base_color_reverse[base], end='')
        print()

    def get_window_label(self, window, positional_candidates, pos_vcf):
        window_label = []
        for pos in range(window[0], window[1]):
            if pos in positional_candidates:
                gt = self.get_combined_gt(self.get_label_of_allele(positional_candidates[pos], pos_vcf))
                window_label.append(gt)
            else:
                window_label.append(0)

        return window_label

    def generate_pileup(self, reads, windows, positional_candidates, vcf_path, train_mode):
        ref_start = max(0, windows[0][0] - CandidateFinderOptions.SAFE_BASES)
        ref_end = windows[-1][1] + CandidateFinderOptions.SAFE_BASES
        reference_sequence = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                       ref_start,
                                                                       ref_end)

        # image generator object
        image_generator = FRIDAY.ImageGenerator(reference_sequence,
                                                self.chromosome_name,
                                                ref_start,
                                                ref_end,
                                                positional_candidates)
        if train_mode:
            vcf_handler = FRIDAY.VCF_handler(vcf_path)
            positional_vcf = vcf_handler.get_positional_vcf_records(self.chromosome_name, ref_start - 20, ref_end + 20)
            image_generator.set_positional_vcf(positional_vcf)

        pileup_images = image_generator.create_window_pileups(windows, reads, train_mode)

        # labels = []
        # for window in windows:
        #     labels.append(self.get_window_label(window, positional_candidates, positional_vcf))
        # print("GOT READS")
        #
        # for i, pileup_image in enumerate(pileup_images):
        #     print(pileup_image.chromosome_name, pileup_image.start_pos, pileup_image.end_pos)
        #     for label in pileup_image.label:
        #         print(label, end='')
        #     print()
        #     for image_row in pileup_image.image:
        #         self.decode_image_row(image_row)
        #
        # exit()
        return pileup_images
