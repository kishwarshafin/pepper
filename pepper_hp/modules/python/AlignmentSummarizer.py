from pepper_hp.build import PEPPER_HP
import itertools
import time
import numpy as np
from operator import itemgetter
from pepper_hp.modules.python.Options import ImageSizeOptions, AlingerOptions, ReadFilterOptions


class AlignmentSummarizer:
    def __init__(self, bam_handler, fasta_handler, chromosome_name, region_start, region_end):
        self.bam_handler = bam_handler
        self.fasta_handler = fasta_handler
        self.chromosome_name = chromosome_name
        self.region_start_position = region_start
        self.region_end_position = region_end

    @staticmethod
    def chunk_images(summary, chunk_size, chunk_overlap):
        chunk_start = 0
        chunk_id = 0
        chunk_end = min(len(summary.genomic_pos), chunk_size)
        images_hp1 = []
        images_hp2 = []
        labels_hp1 = []
        labels_hp2 = []
        labels = []
        positions = []
        chunk_ids = []
        ref_images = []

        while True:
            image_chunk_hp1 = summary.image_hp1[chunk_start:chunk_end]
            image_chunk_hp2 = summary.image_hp2[chunk_start:chunk_end]
            pos_chunk = summary.genomic_pos[chunk_start:chunk_end]
            ref_chunk = summary.ref_image[chunk_start:chunk_end]
            label_chunk_hp1 = [0] * (chunk_end - chunk_start)
            label_chunk_hp2 = [0] * (chunk_end - chunk_start)

            assert (len(image_chunk_hp1) == len(image_chunk_hp2) == len(pos_chunk) == len(label_chunk_hp1) == len(label_chunk_hp2))
            # print(len(image_chunk), len(pos_chunk), len(label_chunk))

            padding_required = chunk_size - len(image_chunk_hp1)
            if padding_required > 0:
                label_chunk_hp1 = label_chunk_hp1 + [0] * padding_required
                label_chunk_hp2 = label_chunk_hp2 + [0] * padding_required
                pos_chunk = pos_chunk + [(-1, -1)] * padding_required
                ref_chunk = ref_chunk + [0] * padding_required
                image_chunk_hp1 = image_chunk_hp1 + [[0.0] * ImageSizeOptions.IMAGE_HEIGHT] * padding_required
                image_chunk_hp2 = image_chunk_hp2 + [[0.0] * ImageSizeOptions.IMAGE_HEIGHT] * padding_required

            assert (len(image_chunk_hp1) == len(image_chunk_hp2) == len(pos_chunk) == len(label_chunk_hp1) == len(label_chunk_hp2) == ImageSizeOptions.SEQ_LENGTH)

            images_hp1.append(image_chunk_hp1)
            images_hp2.append(image_chunk_hp2)
            labels_hp1.append(label_chunk_hp1)
            labels_hp2.append(label_chunk_hp2)
            positions.append(pos_chunk)
            chunk_ids.append(chunk_id)
            ref_images.append(ref_chunk)
            chunk_id += 1

            if chunk_end == len(summary.genomic_pos):
                break

            chunk_start = chunk_end - chunk_overlap
            chunk_end = min(len(summary.genomic_pos), chunk_start + chunk_size)

        return images_hp1, images_hp2, labels_hp1, labels_hp2, positions, chunk_ids, ref_images

    @staticmethod
    def chunk_images_train(summary, chunk_size, chunk_overlap):
        images_hp1 = []
        images_hp2 = []
        labels_hp1 = []
        labels_hp2 = []
        positions = []
        chunk_ids = []
        ref_images = []

        bad_indices = summary.bad_label_positions
        chunk_start = 0
        chunk_id = 0

        for i in range(len(bad_indices)):
            chunk_end = min(chunk_start + chunk_size, bad_indices[i])

            while True:
                if chunk_end - chunk_start != chunk_size:
                    padding_required = chunk_size - (chunk_end - chunk_start)
                    chunk_start -= padding_required
                    if chunk_start < 0:
                        break
                    if i > 0 and chunk_start < bad_indices[i-1]:
                        break
                images_hp1_chunk = summary.image_hp1[chunk_start:chunk_end]
                images_hp2_chunk = summary.image_hp2[chunk_start:chunk_end]
                pos_chunk = summary.genomic_pos[chunk_start:chunk_end]
                ref_chunk = summary.ref_image[chunk_start:chunk_end]
                label_chunk_hp1 = summary.labels_hp1[chunk_start:chunk_end]
                label_chunk_hp2 = summary.labels_hp2[chunk_start:chunk_end]

                assert (len(images_hp1_chunk) == len(images_hp2_chunk) == len(pos_chunk) == len(label_chunk_hp1) == len(label_chunk_hp2) == chunk_size)

                images_hp1.append(images_hp1_chunk)
                images_hp2.append(images_hp2_chunk)
                labels_hp1.append(label_chunk_hp1)
                labels_hp2.append(label_chunk_hp2)
                positions.append(pos_chunk)
                chunk_ids.append(chunk_id)
                ref_images.append(ref_chunk)
                chunk_id += 1

                if chunk_end == bad_indices[i]:
                    break

                chunk_start = chunk_end - chunk_overlap
                chunk_end = min(bad_indices[i], chunk_start + chunk_size)

            chunk_start = chunk_end + 1

        assert(len(images_hp1) == len(images_hp2) == len(labels_hp1) == len(labels_hp2) == len(positions) == len(chunk_ids))

        return images_hp1, images_hp2, labels_hp1, labels_hp2, positions, chunk_ids, ref_images

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    @staticmethod
    def get_overlap_between_ranges(range_a, range_b):
        if range_a[1] > range_b[0]:
            return range_b[0], range_a[1]
        else:
            return None

    def remove_conflicting_regions(self, regions, min_length=ImageSizeOptions.MIN_SEQUENCE_LENGTH,
                                   length_ratio=2.0, overlap_fraction=0.5):
        # reused from medaka's filter_alignments method.
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
    def intersect_intervals(interval_a, intervals):
        # First interval
        l = interval_a[0]
        r = interval_a[1]
        intersection_found = False

        for i in range(0, len(intervals)):
            if intervals[i][0] > r or intervals[i][1] < l:
                continue
            else:
                l = max(l, interval_a[0])
                r = min(r, interval_a[1])
                intersection_found = True

        return l, r, intersection_found

    def reads_to_reference_realignment(self, region_start, region_end, reads):
        # PERFORMS LOCAL REALIGNMENT OF READS TO THE REFERENCE
        if not reads:
            return []

        ref_start = region_start
        ref_end = region_end + AlingerOptions.ALIGNMENT_SAFE_BASES

        ref_sequence = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                 ref_start,
                                                                 ref_end)

        aligner = PEPPER_HP.ReadAligner(ref_start, ref_end, ref_sequence)

        realigned_reads = aligner.align_reads_to_reference(reads)

        # generate_pileup_from_reads.pileup_from_reads(ref_sequence, ref_start, ref_end, realigned_reads)

        return realigned_reads

    @staticmethod
    def range_intersection(intervals):
        left = intervals[0][0]
        right = intervals[0][1]

        found_overlap = False
        for i in range(1, len(intervals)):
            if intervals[i][0] > right or intervals[i][1] < left:
                continue
            else:
                found_overlap = True
                left = max(left, intervals[i][0])
                right = min(right, intervals[i][1])

        if found_overlap:
            return [left, right]
        else:
            return None

    def create_summary(self, truth_bam_h1, truth_bam_h2, train_mode, realignment_flag):
        log_prefix = "[" + self.chromosome_name + ":" + str(self.region_start_position) + "-" \
                     + str(self.region_end_position) + "]"
        all_images_hp1 = []
        all_images_hp2 = []
        all_labels_hp1 = []
        all_labels_hp2 = []
        all_positions = []
        all_image_chunk_ids = []
        all_ref_seq = []

        if train_mode:
            # get the reads from the bam file
            truth_bam_handler_h1 = PEPPER_HP.BAM_handler(truth_bam_h1)
            truth_bam_handler_h2 = PEPPER_HP.BAM_handler(truth_bam_h2)

            truth_reads_h1 = truth_bam_handler_h1.get_reads(self.chromosome_name,
                                                            self.region_start_position,
                                                            self.region_end_position,
                                                            ReadFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                            ReadFilterOptions.MIN_MAPQ,
                                                            ReadFilterOptions.MIN_BASEQ)

            truth_reads_h2 = truth_bam_handler_h2.get_reads(self.chromosome_name,
                                                            self.region_start_position,
                                                            self.region_end_position,
                                                            ReadFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                            ReadFilterOptions.MIN_MAPQ,
                                                            ReadFilterOptions.MIN_BASEQ)

            # do a local realignment of truth reads to reference
            if realignment_flag:
                truth_reads_h1 = self.reads_to_reference_realignment(self.region_start_position,
                                                                     self.region_end_position,
                                                                     truth_reads_h1)
                truth_reads_h2 = self.reads_to_reference_realignment(self.region_start_position,
                                                                     self.region_end_position,
                                                                     truth_reads_h2)

            truth_regions_h1 = []
            for read in truth_reads_h1:
                # start, end, read, is_kept
                truth_regions_h1.append([read.pos, read.pos_end - 1, read,  True])

            truth_regions_h2 = []
            for read in truth_reads_h2:
                # start, end, read, is_kept
                truth_regions_h2.append([read.pos, read.pos_end - 1, read,  True])

            # these are all the regions we will use to generate summaries from.
            # It's important to notice that we need to realign the reads to the reference before we do that.
            truth_regions_h1 = self.remove_conflicting_regions(truth_regions_h1)
            truth_regions_h2 = self.remove_conflicting_regions(truth_regions_h2)
            regions_h1 = []
            for start, end_pos, read, is_kept in truth_regions_h1:
                if is_kept:
                    regions_h1.append([start, end_pos])
            regions_h2 = []
            for start_pos, end_pos, read, is_kept in truth_regions_h2:
                if is_kept:
                    regions_h2.append([start_pos, end_pos])

            truth_regions = []
            for reg_h1 in regions_h1:
                reg = AlignmentSummarizer.range_intersection([reg_h1] + regions_h2)
                if reg is not None:
                    truth_regions.append(reg)

            if not truth_regions:
                # sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " NO TRAINING REGION FOUND.\n"
                #                  + TextColor.END)
                return [], [], [], [], [], [], []

            for region_start, region_end in truth_regions:
                truth_reads_h1 = truth_bam_handler_h1.get_reads(self.chromosome_name,
                                                                self.region_start_position,
                                                                self.region_end_position,
                                                                ReadFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                                ReadFilterOptions.MIN_MAPQ,
                                                                ReadFilterOptions.MIN_BASEQ)

                truth_reads_h2 = truth_bam_handler_h2.get_reads(self.chromosome_name,
                                                                self.region_start_position,
                                                                self.region_end_position,
                                                                ReadFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                                ReadFilterOptions.MIN_MAPQ,
                                                                ReadFilterOptions.MIN_BASEQ)
                if len(truth_reads_h1) > 1 or len(truth_reads_h2) > 1:
                    continue

                ref_start = region_start
                ref_end = region_end + 1
                # ref_seq should contain region_end_position base
                ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                    ref_start,
                                                                    ref_end)

                read_start = max(0, region_start)
                read_end = region_end
                all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                       read_start,
                                                       read_end,
                                                       ReadFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                       ReadFilterOptions.MIN_MAPQ,
                                                       ReadFilterOptions.MIN_BASEQ)
                total_reads = len(all_reads)

                if total_reads == 0:
                    continue

                if total_reads > AlingerOptions.MAX_READS_IN_REGION:
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

                # sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " TOTAL " + str(total_reads)
                #                  + " READS FOUND.\n" + TextColor.END)

                start_time = time.time()

                if realignment_flag:
                    all_reads = self.reads_to_reference_realignment(read_start,
                                                                    read_end,
                                                                    all_reads)
                    # sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " REALIGNMENT OF TOTAL "
                    #                  + str(total_reads) + " READS TOOK: " + str(round(time.time()-start_time, 5))
                    #                  + " secs\n" + TextColor.END)

                summary_generator = PEPPER_HP.SummaryGenerator(ref_seq,
                                                               self.chromosome_name,
                                                               ref_start,
                                                               ref_end)

                summary_generator.generate_train_summary(all_reads,
                                                         region_start,
                                                         region_end,
                                                         truth_reads_h1[0],
                                                         truth_reads_h2[0])

                images_hp1, images_hp2, labels_hp1, labels_hp2, positions, chunk_ids, ref_seqs = \
                    self.chunk_images_train(summary_generator,
                                            chunk_size=ImageSizeOptions.SEQ_LENGTH,
                                            chunk_overlap=ImageSizeOptions.SEQ_OVERLAP)

                all_images_hp1.extend(images_hp1)
                all_images_hp2.extend(images_hp2)
                all_labels_hp1.extend(labels_hp1)
                all_labels_hp2.extend(labels_hp2)
                all_positions.extend(positions)
                all_image_chunk_ids.extend(chunk_ids)
                all_ref_seq.extend(ref_seqs)
        else:
            # HERE REALIGN THE READS TO THE REFERENCE THEN GENERATE THE SUMMARY TO GET A POLISHED HAPLOTYPE
            read_start = max(0, self.region_start_position)
            read_end = self.region_end_position

            all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                   read_start,
                                                   read_end,
                                                   ReadFilterOptions.INCLUDE_SUPPLEMENTARY,
                                                   ReadFilterOptions.MIN_MAPQ,
                                                   ReadFilterOptions.MIN_BASEQ)

            total_reads = len(all_reads)

            if total_reads == 0:
                return [], [], [], [], [], [], []

            if total_reads > AlingerOptions.MAX_READS_IN_REGION:
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

            # sys.stderr.write(TextColor.PURPLE + "INFO: " + log_prefix + " TOTAL " + str(total_reads) + " READS FOUND\n"
            #                  + TextColor.END)

            if realignment_flag:
                start_time = time.time()
                all_reads = self.reads_to_reference_realignment(self.region_start_position,
                                                                self.region_end_position,
                                                                all_reads)
                # sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " REALIGNMENT OF TOTAL " + str(total_reads)
                #                 + " READS TOOK: " + str(round(time.time()-start_time, 5)) + " secs\n" + TextColor.END)

            # ref_seq should contain region_end_position base
            ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                self.region_start_position,
                                                                self.region_end_position + 1)

            summary_generator = PEPPER_HP.SummaryGenerator(ref_seq,
                                                           self.chromosome_name,
                                                           self.region_start_position,
                                                           self.region_end_position)

            summary_generator.generate_summary(all_reads,
                                               self.region_start_position,
                                               self.region_end_position)


            # images_hp1, images_hp2, labels, positions, chunk_ids, ref_images
            images_hp1, images_hp2, labels_hp1, labels_hp2, positions, chunk_ids, ref_seqs = \
                self.chunk_images(summary_generator,
                                  chunk_size=ImageSizeOptions.SEQ_LENGTH,
                                  chunk_overlap=ImageSizeOptions.SEQ_OVERLAP)

            all_images_hp1.extend(images_hp1)
            all_images_hp2.extend(images_hp2)
            all_labels_hp1.extend(labels_hp1)
            all_labels_hp2.extend(labels_hp2)
            all_positions.extend(positions)
            all_image_chunk_ids.extend(chunk_ids)
            all_ref_seq.extend(ref_seqs)

        assert(len(all_images_hp1) == len(all_images_hp2) == len(all_labels_hp1) == len(all_labels_hp2) == len(all_image_chunk_ids) == len(all_ref_seq))

        return all_images_hp1, all_images_hp2, all_labels_hp1, all_labels_hp2, all_positions, all_image_chunk_ids, all_ref_seq
