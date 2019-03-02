from build import HELEN
import numpy as np
import math
from modules.python.ActiveRegionFinder import ActiveRegionFinder, ActiveRegionOptions
from modules.python.DeBruijnHaplotyper import DeBruijnHaplotyper
from modules.python.Options import AlingerOptions, CandidateFinderOptions, DeBruijnGraphOptions


class RegionBasedHaplotypes:
    def __init__(self, haplotypes, region_start, region_end):
        self.haplotypes = haplotypes
        self.region_start = region_start
        self.region_end = region_end
        self.min_read_start = None
        self.max_read_end = None
        self.reads = []

    def assign_read(self, read):
        if self.min_read_start is None or self.max_read_end is None:
            self.min_read_start = read.pos
            self.max_read_end = read.pos_end

        self.min_read_start = min(self.min_read_start, read.pos)
        self.max_read_end = max(self.max_read_end, read.pos_end)
        self.reads.append(read)

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))


class LocalAssembler:
    def __init__(self, bam_handler, fasta_handler, chromosome_name, region_start, region_end):
        self.bam_handler = bam_handler
        self.fasta_handler = fasta_handler
        self.chromosome_name = chromosome_name
        self.region_start_position = region_start
        self.region_end_position = region_end

    def perform_local_alignment(self, region_with_reads):
        if not region_with_reads.reads:
            return []
        ref_start = min(region_with_reads.min_read_start, region_with_reads.region_start) - AlingerOptions.ALIGNMENT_SAFE_BASES
        ref_end = max(region_with_reads.max_read_end, region_with_reads.region_end) + AlingerOptions.ALIGNMENT_SAFE_BASES

        # ref_start = region_with_reads.region_start - AlingerOptions.ALIGNMENT_SAFE_BASES
        # ref_end = region_with_reads.region_end + AlingerOptions.ALIGNMENT_SAFE_BASES

        if ref_end <= region_with_reads.region_end:
            return region_with_reads.reads
        else:
            ref_suffix = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                   region_with_reads.region_end,
                                                                   ref_end)

        ref_prefix = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                               ref_start,
                                                               region_with_reads.region_start)
        ref = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                        region_with_reads.region_start,
                                                        region_with_reads.region_end)
        ref_seq = ref_prefix + ref + ref_suffix
        haplotypes = [ref_prefix + hap + ref_suffix for hap in region_with_reads.haplotypes]

        aligner = FRIDAY.ReadAligner(ref_start, ref_end, ref_seq)

        haplotypes = sorted(set(haplotypes))

        if not haplotypes or haplotypes == [ref_seq]:
            return region_with_reads.reads

        realigned_reads = aligner.align_reads(haplotypes, region_with_reads.reads)

        return realigned_reads

    def perform_local_assembly(self, downsample_rate, perform_alignment=True):
        # get the reads from the bam file
        all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                               self.region_start_position,
                                               self.region_end_position,
                                               CandidateFinderOptions.MIN_MAP_QUALITY,
                                               DeBruijnGraphOptions.MIN_BASE_QUALITY)

        reads_to_keep = min(int(math.ceil(len(all_reads) * downsample_rate)), int(AlingerOptions.MAX_READS_IN_REGION))

        if len(all_reads) > reads_to_keep:
            # https://github.com/google/nucleus/blob/master/nucleus/util/utils.py
            # reservoir_sample method utilized here
            random = np.random.RandomState(AlingerOptions.RANDOM_SEED)
            sample = []
            for i, read in enumerate(all_reads):
                if len(sample) < AlingerOptions.MAX_READS_IN_REGION:
                    sample.append(read)
                else:
                    j = random.randint(0, i + 1)
                    if j < AlingerOptions.MAX_READS_IN_REGION:
                        sample[j] = read

            all_reads = sample

        if perform_alignment is False:
            return all_reads

        # find active regions
        active_region_finder = ActiveRegionFinder(self.chromosome_name,
                                                  self.region_start_position,
                                                  self.region_end_position,
                                                  self.fasta_handler)

        # mapq threshold is checked in this
        active_regions = active_region_finder.find_active_region(all_reads)
        assembly_active_regions = []
        possible_regions = []

        for active_region in active_regions:
            start_pos, end_pos = active_region

            if end_pos - start_pos > ActiveRegionOptions.MAX_REGION_SIZE:
                continue

            db_graph = DeBruijnHaplotyper(self.fasta_handler,
                                          self.chromosome_name,
                                          start_pos,
                                          end_pos)

            reference_sequence, haplotypes = db_graph.find_haplotypes(all_reads)

            # print(active_region)
            # for hap in set(haplotypes):
            #     print(hap)
            if haplotypes:
                assembly_active_regions.append(RegionBasedHaplotypes(haplotypes, start_pos, end_pos))
                possible_regions.append((start_pos, end_pos))

        if not possible_regions:
            return all_reads
        # now we have the list that we filtered at the beginning of this script
        realigned_reads = list()
        for read in all_reads:
            read_range = (read.pos, read.pos_end)
            overlapping_lengths = [RegionBasedHaplotypes.overlap_length_between_ranges(region, read_range)
                                   for region in possible_regions]
            max_length = max(overlapping_lengths)
            if max_length <= 0:
                realigned_reads.append(read)
                continue

            max_window_index = max(range(len(possible_regions)), key=lambda i: overlapping_lengths[i])
            assembly_active_regions[max_window_index].assign_read(read)

        for active_region in assembly_active_regions:
            realigned_reads.extend(self.perform_local_alignment(active_region))

        realigned_reads = sorted(realigned_reads)

        return realigned_reads
