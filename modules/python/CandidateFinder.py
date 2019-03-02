from build import HELEN
from modules.python.Options import CandidateFinderOptions, ImageSizeOptions


class CandidateFinder:
    def __init__(self, fasta_handler, contig, start, end):
        self.fasta_handler = fasta_handler
        self.contig = contig
        self.region_start = start
        self.region_end = end

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    def find_candidates(self, reads):
        ref_start = max(0, self.region_start - (CandidateFinderOptions.SAFE_BASES * 2))
        ref_end = self.region_end + (CandidateFinderOptions.SAFE_BASES * 2)

        reference_sequence = self.fasta_handler.get_reference_sequence(self.contig,
                                                                       ref_start,
                                                                       ref_end)
        # candidate finder objects
        candidate_finder = FRIDAY.CandidateFinder(reference_sequence,
                                                  self.contig,
                                                  max(0, self.region_start - CandidateFinderOptions.SAFE_BASES),
                                                  self.region_end + CandidateFinderOptions.SAFE_BASES,
                                                  ref_start,
                                                  ref_end)

        # find candidates
        candidate_positions, candidate_map = candidate_finder.find_candidates(reads)

        filtered_positions = []
        for candidate_position in candidate_positions:
            if self.region_start <= candidate_position <= self.region_end:
                filtered_positions.append(candidate_position)

        return filtered_positions, candidate_map

    @staticmethod
    def get_windows_from_candidates(candidate_positions):
        current_window_st, current_window_end = -1, -1
        all_windows = []
        for candidate_position in sorted(candidate_positions):
            if current_window_st == -1:
                current_window_st = candidate_position - ImageSizeOptions.MIN_BASES_ON_LEFT
                current_window_end = candidate_position + ImageSizeOptions.BASES_ON_RIGHT
                all_windows.append((current_window_st, current_window_end))
            elif candidate_position < current_window_end - ImageSizeOptions.MIN_BASES_ON_LEFT:
                continue
            else:
                current_window_st = candidate_position - ImageSizeOptions.MIN_BASES_ON_LEFT
                current_window_end = candidate_position + ImageSizeOptions.BASES_ON_RIGHT
                all_windows.append((current_window_st, current_window_end))

        return sorted(all_windows)

