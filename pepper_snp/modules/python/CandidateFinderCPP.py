from collections import defaultdict
from pepper_snp.build import PEPPER_SNP
from pepper_snp.modules.python.Options import ImageSizeOptions, ReadFilterOptions, CandidateFinderOptions


class CandidateFinderCPP:
    def __init__(self, contig, start, end):
        self.contig = contig
        self.region_start = start
        self.region_end = end

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    def find_candidates(self, bam_file_path, fasta_file_path, contig_name, region_start, region_end, all_positions, all_indicies, all_predictions, set_profile):
        bam_handler = PEPPER_SNP.BAM_handler(bam_file_path)
        fasta_handler = PEPPER_SNP.FASTA_handler(fasta_file_path)
        all_reads = bam_handler.get_reads(contig_name,
                                          region_start,
                                          region_end,
                                          ReadFilterOptions.INCLUDE_SUPPLEMENTARY,
                                          ReadFilterOptions.MIN_MAPQ,
                                          ReadFilterOptions.MIN_BASEQ)

        ref_start = max(0, self.region_start - (CandidateFinderOptions.SAFE_BASES * 2))
        ref_end = self.region_end + (CandidateFinderOptions.SAFE_BASES * 2)

        reference_sequence = fasta_handler.get_reference_sequence(self.contig,
                                                                  ref_start,
                                                                  ref_end).upper()

        # candidate finder objects
        candidate_finder = PEPPER_SNP.CandidateFinder(reference_sequence,
                                                      contig_name,
                                                      max(0, region_start - CandidateFinderOptions.SAFE_BASES),
                                                      region_end + CandidateFinderOptions.SAFE_BASES,
                                                      ref_start,
                                                      ref_end)
        # find candidates
        positional_candidate_list = candidate_finder.find_candidates(all_reads, all_positions, all_indicies, all_predictions, set_profile)
        positional_candidate_list = sorted(positional_candidate_list)
        # print(positional_candidate_list)

        positional_candidate_map = defaultdict(list)
        for positional_candidate in positional_candidate_list:
            for candidate in positional_candidate.candidates:
                if region_start <= candidate.pos_start and candidate.pos_end <= region_end:
                    positional_candidate_map[positional_candidate.pos_start].append(candidate)

        return positional_candidate_map
