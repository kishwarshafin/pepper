from pepper_variant.build import PEPPER_VARIANT
import itertools
import numpy as np
from operator import itemgetter
from pepper_variant.modules.python.Options import ImageSizeOptions, AlingerOptions, ReadFilterOptions, TruthFilterOptions


class AlignmentSummarizer:
    """

    """
    def __init__(self, bam_handler, fasta_handler, chromosome_name, region_start, region_end):
        self.bam_handler = bam_handler
        self.fasta_handler = fasta_handler
        self.chromosome_name = chromosome_name
        self.region_start_position = region_start
        self.region_end_position = region_end

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    @staticmethod
    def get_overlap_between_ranges(range_a, range_b):
        if range_a[1] > range_b[0]:
            return range_b[0], range_a[1]
        else:
            return None

    def remove_conflicting_regions(self, regions, min_length=ImageSizeOptions.MIN_SEQUENCE_LENGTH, length_ratio=2.0, overlap_fraction=0.5):
        for reg_a, reg_b in itertools.combinations(regions, 2):
            el1, el2 = sorted((reg_a, reg_b), key=itemgetter(0))
            overlap = self.get_overlap_between_ranges(el1, el2)

            if overlap is None:
                continue
            ovlp_start, ovlp_end = overlap
            s, l = sorted((reg_a, reg_b), key=lambda element: (element[1] - element[0]))

            length_ratio_ij = (l[1] - l[0]) / max(1, (s[1] - s[0]))
            overlap_fraction_ij = (ovlp_end - ovlp_start) / max(1, (s[1] - s[0]))
            # 4 cases
            if length_ratio_ij < length_ratio:  # we don't trust one more than the other
                if overlap_fraction_ij >= overlap_fraction:
                    # 1) they overlap a lot; we have significant ambiguity, discard both
                    s[3] = False
                    l[3] = False
                else:
                    # 2) they overlap a little; just trim overlapping portions of both alignments
                    el1[1] = ovlp_start
                    el2[0] = ovlp_end
            else:  # we trust the longer one more than the shorter one
                if overlap_fraction_ij >= overlap_fraction:
                    # 3) they overlap a lot; discard shorter alignment
                    s[3] = False
                else:
                    # 4) they overlap a little; trim overlapping portion of shorter alignment
                    el2[0] = ovlp_end

        # trim starts and ends if needed
        for al in regions:
            al[0] = max(self.region_start_position, al[0])
            al[1] = min(self.region_end_position, al[1])
        # do filtering
        filtered_alignments = [al for al in regions
                               if (al[3] and al[1] - al[0] >= min_length)]
        filtered_alignments.sort(key=itemgetter(0))

        return filtered_alignments

    @staticmethod
    def range_intersection(intervals):
        left = intervals[0][0]
        right = intervals[0][1]
        read_hp1 = intervals[0][2]
        read_hp2 = None

        longest_interval = 0
        found_overlap = False
        for i in range(1, len(intervals)):
            if intervals[i][0] > right or intervals[i][1] < left:
                continue
            else:
                found_overlap = True
                overlap_length = intervals[i][1] - intervals[i][0]
                if overlap_length > longest_interval:
                    longest_interval = overlap_length
                    left = intervals[i][0]
                    right = intervals[i][1]
                    read_hp2 = intervals[i][2]

        if found_overlap:
            return [left, right, read_hp1, read_hp2]
        else:
            return None

    @staticmethod
    def range_intersection_bed(interval, bed_intervals):
        left = interval[0]
        right = interval[1]
        read_hp1 = interval[2]
        read_hp2 = interval[3]

        intervals = []
        for i in range(0, len(bed_intervals)):
            bed_left = bed_intervals[i][0]
            bed_right = bed_intervals[i][1]

            if bed_right < left:
                continue
            elif bed_left > right:
                continue
            else:
                left_bed = max(left, bed_left)
                right_bed = min(right, bed_right)
                intervals.append([left_bed, right_bed, read_hp1, read_hp2])

        return intervals

    def create_summary(self, truth_bam_handler_h1, truth_bam_handler_h2, train_mode, downsample_rate, bed_list):
        all_images = []
        all_labels = []
        all_positions = []
        all_index = []
        all_chunk_ids = []

        if train_mode:
            # get the reads from the bam file
            truth_reads_h1 = truth_bam_handler_h1.get_reads(self.chromosome_name,
                                                            self.region_start_position,
                                                            self.region_end_position,
                                                            TruthFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                            TruthFilterOptions.MIN_MAPQ,
                                                            TruthFilterOptions.MIN_BASEQ)

            truth_reads_h2 = truth_bam_handler_h2.get_reads(self.chromosome_name,
                                                            self.region_start_position,
                                                            self.region_end_position,
                                                            TruthFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                            TruthFilterOptions.MIN_MAPQ,
                                                            TruthFilterOptions.MIN_BASEQ)

            truth_regions_h1 = []
            for read in truth_reads_h1:
                # start, end, read, is_kept
                truth_regions_h1.append([read.pos, read.pos_end - 1, read,  True])

            truth_regions_h2 = []
            for read in truth_reads_h2:
                # start, end, read, is_kept
                truth_regions_h2.append([read.pos, read.pos_end - 1, read,  True])

            # these are all the regions we will use to generate summaries from.
            truth_regions_h1 = self.remove_conflicting_regions(truth_regions_h1)
            truth_regions_h2 = self.remove_conflicting_regions(truth_regions_h2)

            regions_h1 = []
            for start, end_pos, read, is_kept in truth_regions_h1:
                if is_kept:
                    regions_h1.append([start, end_pos, read])

            regions_h2 = []
            for start_pos, end_pos, read, is_kept in truth_regions_h2:
                if is_kept:
                    regions_h2.append([start_pos, end_pos, read])

            truth_regions = []
            for reg_h1 in regions_h1:
                reg = AlignmentSummarizer.range_intersection([reg_h1] + regions_h2)
                if reg is not None:
                    truth_regions.append(reg)

            # now intersect with bed file
            if bed_list is not None:
                intersected_truth_regions = []
                for region in truth_regions:
                    if self.chromosome_name in bed_list.keys():
                        reg = AlignmentSummarizer.range_intersection_bed(region, bed_list[self.chromosome_name])
                        intersected_truth_regions.extend(reg)
                truth_regions = intersected_truth_regions

            if not truth_regions:
                # sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " NO TRAINING REGION FOUND.\n"
                #                  + TextColor.END)
                return None

            chunk_id_start = 0

            for region_start, region_end, truth_read_h1, truth_read_h2 in truth_regions:
                if not truth_reads_h1 or not truth_reads_h2:
                    continue

                ref_start = max(0, region_start)
                ref_end = region_end

                # ref_seq should contain region_end_position base
                ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                    ref_start,
                                                                    ref_end + 1)

                read_start = max(0, region_start)
                read_end = region_end
                all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                       read_start,
                                                       read_end,
                                                       ReadFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                       ReadFilterOptions.MIN_MAPQ,
                                                       ReadFilterOptions.MIN_BASEQ)

                total_reads = len(all_reads)
                total_allowed_reads = int(min(AlingerOptions.MAX_READS_IN_REGION, downsample_rate * total_reads))

                if total_reads > total_allowed_reads:
                    # https://github.com/google/nucleus/blob/master/nucleus/util/utils.py
                    # reservoir_sample method utilized here
                    random = np.random.RandomState(AlingerOptions.RANDOM_SEED)
                    sample = []
                    for i, read in enumerate(all_reads):
                        if len(sample) < total_allowed_reads:
                            sample.append(read)
                        else:
                            j = random.randint(0, i + 1)
                            if j < total_allowed_reads:
                                sample[j] = read
                    all_reads = sample

                total_reads = len(all_reads)

                if total_reads == 0:
                    continue

                regional_summary = PEPPER_VARIANT.RegionalSummaryGenerator(region_start, region_end, ref_seq)

                regional_summary.generate_max_insert_summary(all_reads)

                regional_summary.generate_labels(truth_reads_h1[0],
                                                 truth_reads_h2[0])

                regional_image_summary = regional_summary.generate_summary(all_reads,
                                                                           ImageSizeOptions.SEQ_OVERLAP,
                                                                           ImageSizeOptions.SEQ_LENGTH,
                                                                           ImageSizeOptions.IMAGE_HEIGHT,
                                                                           chunk_id_start,
                                                                           train_mode)




                #############################
                all_images.extend(regional_image_summary.chunked_image_matrix)
                all_labels.extend(regional_image_summary.chunked_labels)
                all_positions.extend(regional_image_summary.chunked_positions)
                all_index.extend(regional_image_summary.chunked_index)
                all_chunk_ids.extend(regional_image_summary.chunked_ids)

                if len(regional_image_summary.chunked_ids) > 0:
                    chunk_id_start = max(regional_image_summary.chunked_ids) + 1

                # summary_generator = PEPPER_VARIANT.SummaryGenerator(ref_seq,
                #                                                     self.chromosome_name,
                #                                                     ref_start,
                #                                                     ref_end)
                #
                # summary_generator.generate_train_summary(all_reads,
                #                                          region_start,
                #                                          region_end,
                #                                          truth_reads_h1[0],
                #                                          truth_reads_h2[0])
                #
                # image_summary = summary_generator.chunk_image_train(ImageSizeOptions.SEQ_LENGTH,
                #                                                     ImageSizeOptions.SEQ_OVERLAP,
                #                                                     ImageSizeOptions.IMAGE_HEIGHT,
                #                                                     chunk_id_start)
                #
                # for i in range(0, len(image_summary.positions)):
                #     for j in range(0, len(image_summary.positions[i])):
                #         if image_summary.positions[i][j][0] != regional_image_summary.chunked_positions[i][j]:
                #             print("NOT EQUAL", i, j, image_summary.positions[i][j][0], regional_image_summary.chunked_positions[i][j])
                # print("PASSED POSITIONS")
                #
                # for i in range(0, len(image_summary.positions)):
                #     for j in range(0, len(image_summary.positions[i])):
                #         if image_summary.positions[i][j][1] != regional_image_summary.chunked_index[i][j]:
                #             print("NOT EQUAL", i, j, image_summary.positions[i][j][1], regional_image_summary.chunked_index[i][j])
                # print("PASSED INDEX")
                #
                # for i in range(0, len(image_summary.chunk_ids)):
                #     if image_summary.chunk_ids[i] != regional_image_summary.chunked_ids[i]:
                #         print("NOT EQUAL", i, image_summary.chunk_ids[i], regional_image_summary.chunked_ids[i])
                # print("PASSED CHUNK IDs")
                #
                # for i in range(0, len(image_summary.labels)):
                #     for j in range(0, len(image_summary.labels[i])):
                #         if image_summary.labels[i][j] != regional_image_summary.chunked_labels[i][j]:
                #             print("NOT EQUAL", i, j, image_summary.labels[i][j][1], regional_image_summary.chunked_labels[i][j])
                # print("PASSED LABELS")

                # for i in range(0, len(image_summary.images)):
                #     for j in range(0, len(image_summary.images[i])):
                #         for k in range(0, len(image_summary.images[i][j])):
                #             if image_summary.images[i][j][k] != regional_image_summary.chunked_image_matrix[i][j][k]:
                #                 print("NOT EQUAL", i, j, image_summary.images[i][j][k], regional_image_summary.chunked_image_matrix[i][j][k])
                # print("PASSED IMAGES")
        else:
            read_start = max(0, self.region_start_position)
            read_end = self.region_end_position

            all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                   read_start,
                                                   read_end,
                                                   ReadFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                   ReadFilterOptions.MIN_MAPQ,
                                                   ReadFilterOptions.MIN_BASEQ)

            total_reads = len(all_reads)
            total_allowed_reads = int(min(AlingerOptions.MAX_READS_IN_REGION, downsample_rate * total_reads))

            if total_reads > total_allowed_reads:
                # sys.stderr.write("INFO: " + log_prefix + "HIGH COVERAGE CHUNK: " + str(total_reads) + " Reads.\n")
                # https://github.com/google/nucleus/blob/master/nucleus/util/utils.py
                # reservoir_sample method utilized here
                random = np.random.RandomState(AlingerOptions.RANDOM_SEED)
                sample = []
                for i, read in enumerate(all_reads):
                    if len(sample) < total_allowed_reads:
                        sample.append(read)
                    else:
                        j = random.randint(0, i + 1)
                        if j < total_allowed_reads:
                            sample[j] = read
                all_reads = sample

            total_reads = len(all_reads)

            if total_reads == 0:
                return [], [], [], [], []

            # ref_seq should contain region_end_position base
            ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                self.region_start_position,
                                                                self.region_end_position + 1)

            chunk_id_start = 0

            regional_summary = PEPPER_VARIANT.RegionalSummaryGenerator(self.region_start_position, self.region_end_position, ref_seq)

            regional_summary.generate_max_insert_summary(all_reads)

            regional_image_summary = regional_summary.generate_summary(all_reads,
                                                                       ImageSizeOptions.SEQ_OVERLAP,
                                                                       ImageSizeOptions.SEQ_LENGTH,
                                                                       ImageSizeOptions.IMAGE_HEIGHT,
                                                                       chunk_id_start,
                                                                       train_mode)

            all_images.extend(regional_image_summary.chunked_image_matrix)
            all_labels.extend(regional_image_summary.chunked_labels)
            all_positions.extend(regional_image_summary.chunked_positions)
            all_index.extend(regional_image_summary.chunked_index)
            all_chunk_ids.extend(regional_image_summary.chunked_ids)



            ############ Test
            summary_generator = PEPPER_VARIANT.SummaryGenerator(ref_seq,
                                                                self.chromosome_name,
                                                                self.region_start_position,
                                                                self.region_end_position)
            summary_generator.generate_summary(all_reads,
                                               self.region_start_position,
                                               self.region_end_position)

            image_summary = summary_generator.chunk_image(ImageSizeOptions.SEQ_LENGTH,
                                                          ImageSizeOptions.SEQ_OVERLAP,
                                                          ImageSizeOptions.IMAGE_HEIGHT)

            for i in range(0, len(image_summary.positions)):
                for j in range(0, len(image_summary.positions[i])):
                    if image_summary.positions[i][j][0] != regional_image_summary.chunked_positions[i][j]:
                        print("NOT EQUAL", i, j, image_summary.positions[i][j][0], regional_image_summary.chunked_positions[i][j])
            print("PASSED POSITIONS")

            for i in range(0, len(image_summary.positions)):
                for j in range(0, len(image_summary.positions[i])):
                    if image_summary.positions[i][j][1] != regional_image_summary.chunked_index[i][j]:
                        print("NOT EQUAL", i, j, image_summary.positions[i][j][1], regional_image_summary.chunked_index[i][j])
            print("PASSED INDEX")

            for i in range(0, len(image_summary.chunk_ids)):
                if image_summary.chunk_ids[i] != regional_image_summary.chunked_ids[i]:
                    print("NOT EQUAL", i, image_summary.chunk_ids[i], regional_image_summary.chunked_ids[i])
            print("PASSED CHUNK IDs")

            for i in range(0, len(image_summary.labels)):
                for j in range(0, len(image_summary.labels[i])):
                    if image_summary.labels[i][j] != regional_image_summary.chunked_labels[i][j]:
                        print("NOT EQUAL", i, j, image_summary.labels[i][j][1], regional_image_summary.chunked_labels[i][j])
            print("PASSED LABELS")

            for i in range(0, len(image_summary.images)):
                for j in range(0, len(image_summary.images[i])):
                    for k in range(0, len(image_summary.images[i][j])):
                        if image_summary.images[i][j][k] != regional_image_summary.chunked_image_matrix[i][j][k]:
                            print("NOT EQUAL", i, j, image_summary.images[i][j][k], regional_image_summary.chunked_image_matrix[i][j][k])
            print("PASSED IMAGES")

        return all_images, all_labels, all_positions, all_index, all_chunk_ids
