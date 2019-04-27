from build import HELEN
from modules.python.Options import AlingerOptions


class Haplotype2Variant:
    def __init__(self, fasta_file_path, contig, region_start, region_end):
        self.contig = contig
        self.region_start = region_start
        self.region_end = region_end
        self.fasta_handler = HELEN.FASTA_handler(fasta_file_path)

    def generate_records_from_haplotypes(self, haplotypes):
        ref_start = self.region_start - AlingerOptions.ALIGNMENT_SAFE_BASES
        ref_end = self.region_end + AlingerOptions.ALIGNMENT_SAFE_BASES
        ref_sequence = self.fasta_handler.get_reference_sequence(self.contig,
                                                                 ref_start,
                                                                 ref_end)

        aligner = HELEN.ReadAligner(ref_start, ref_end, ref_sequence)
        haplotype_reads = aligner.align_haplotypes_to_reference(haplotypes)

        candidate_finder = HELEN.CandidateFinder(ref_sequence,
                                                 self.contig,
                                                 self.region_start,
                                                 self.region_end,
                                                 ref_start,
                                                 ref_end)
        all_called_candidate_variants = candidate_finder.find_candidates(haplotype_reads)

        variant_set = list()
        for variant in all_called_candidate_variants:
            if self.region_start <= variant.pos <= self.region_end:
                genotype = [0, 0]
                if variant.alt2 == '.':
                    alternate_alleles = [variant.alt1]
                    if variant.alt1_count == 1:
                        genotype = [0, 1]
                    elif variant.alt1_count == 2:
                        genotype = [1, 1]
                else:
                    alternate_alleles = [variant.alt1, variant.alt2]
                    if variant.alt1_count == 1 and variant.alt2_count == 1:
                        genotype = [1, 2]

                variant_set.append((variant.chromosome_name, variant.pos, variant.pos_end, variant.ref,
                                    alternate_alleles, genotype))

        return variant_set
