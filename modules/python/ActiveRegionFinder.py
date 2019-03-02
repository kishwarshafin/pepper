from build import HELEN
from modules.python.Options import ActiveRegionOptions


class ActiveRegionFinder:
    def __init__(self, contig, start, end, fasta_handler):
        self.contig = contig
        self.region_start = max(0, start - ActiveRegionOptions.REGION_EXPANSION)
        self.region_end = end + ActiveRegionOptions.REGION_EXPANSION
        self.reference_sequence = fasta_handler.get_reference_sequence(self.contig,
                                                                       self.region_start,
                                                                       self.region_end)

    def find_active_region(self, reads):
        # find the active region
        active_region_finder = FRIDAY.ActiveRegionFinder(self.reference_sequence,
                                                         self.contig,
                                                         self.region_start,
                                                         self.region_end)

        # find active regions
        active_regions = active_region_finder.find_active_region(reads)

        return active_regions
