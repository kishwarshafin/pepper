from pepper.build import PEPPER
import itertools
import time
import sys
import numpy as np
from operator import itemgetter
from pepper.modules.python.Options import ImageSizeOptions, AlingerOptions


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
        images = []
        labels = []
        positions = []
        chunk_ids = []

        while True:
            image_chunk = summary.image[chunk_start:chunk_end]
            pos_chunk = summary.genomic_pos[chunk_start:chunk_end]
            label_chunk = [0] * (chunk_end - chunk_start)

            # assert (len(image_chunk) == len(pos_chunk) == len(label_chunk))
            # print(len(image_chunk), len(pos_chunk), len(label_chunk))

            padding_required = chunk_size - len(image_chunk)
            if padding_required > 0:
                label_chunk = label_chunk + [0] * padding_required
                pos_chunk = pos_chunk + [(-1, -1)] * padding_required
                image_chunk = image_chunk + [[0.0] * ImageSizeOptions.IMAGE_HEIGHT] * padding_required

            # assert (len(image_chunk) == len(pos_chunk) == len(label_chunk) == ImageSizeOptions.SEQ_LENGTH)

            images.append(image_chunk)
            labels.append(label_chunk)
            positions.append(pos_chunk)
            chunk_ids.append(chunk_id)
            chunk_id += 1

            if chunk_end == len(summary.genomic_pos):
                break

            chunk_start = chunk_end - chunk_overlap
            chunk_end = min(len(summary.genomic_pos), chunk_start + chunk_size)

        return images, labels, positions, chunk_ids

    @staticmethod
    def chunk_images_train(summary, chunk_size, chunk_overlap):
        images = []
        labels = []
        positions = []
        chunk_ids = []

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

                image_chunk = summary.image[chunk_start:chunk_end]
                pos_chunk = summary.genomic_pos[chunk_start:chunk_end]
                label_chunk = summary.labels[chunk_start:chunk_end]

                # assert (len(image_chunk) == len(pos_chunk) == len(label_chunk) == chunk_size)

                images.append(image_chunk)
                labels.append(label_chunk)
                positions.append(pos_chunk)
                chunk_ids.append(chunk_id)
                chunk_id += 1

                if chunk_end == bad_indices[i]:
                    break

                chunk_start = chunk_end - chunk_overlap
                chunk_end = min(bad_indices[i], chunk_start + chunk_size)

            chunk_start = chunk_end + 1

        # assert(len(images) == len(labels) == len(positions) == len(chunk_ids))

        return images, labels, positions, chunk_ids

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

    def reads_to_reference_realignment(self, region_start, region_end, reads):
        # PERFORMS LOCAL REALIGNMENT OF READS TO THE REFERENCE
        if not reads:
            return []

        ref_start = region_start
        ref_end = region_end + AlingerOptions.ALIGNMENT_SAFE_BASES

        ref_sequence = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                 ref_start,
                                                                 ref_end)

        aligner = PEPPER.ReadAligner(ref_start, ref_end, ref_sequence)

        realigned_reads = aligner.align_reads_to_reference(reads)

        # generate_pileup_from_reads.pileup_from_reads(ref_sequence, ref_start, ref_end, realigned_reads)

        return realigned_reads

    def create_summary(self, truth_bam_handler, train_mode, downsample_rate, realignment_flag=True):
        log_prefix = "[" + self.chromosome_name + ":" + str(self.region_start_position) + "-" \
                     + str(self.region_end_position) + "]"
        all_images = []
        all_labels = []
        all_positions = []
        all_image_chunk_ids = []

        if train_mode:
            # get the reads from the bam file
            include_supplementary = True
            min_mapq = 60
            min_baseq = 0
            truth_reads = truth_bam_handler.get_reads(self.chromosome_name,
                                                      self.region_start_position,
                                                      self.region_end_position,
                                                      include_supplementary,
                                                      60,
                                                      0)

            # do a local realignment of truth reads to reference
            if realignment_flag:
                truth_reads = self.reads_to_reference_realignment(self.region_start_position,
                                                                  self.region_end_position,
                                                                  truth_reads)

            truth_regions = []
            for read in truth_reads:
                # start, end, read, is_kept, is_h1
                truth_regions.append([read.pos, read.pos_end - 1, read,  True])

            # these are all the regions we will use to generate summaries from.
            # It's important to notice that we need to realign the reads to the reference before we do that.
            truth_regions = self.remove_conflicting_regions(truth_regions)

            if not truth_regions:
                # sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " NO TRAINING REGION FOUND.\n"
                #                  + TextColor.END)
                return [], [], [], []

            for region in truth_regions:
                region_start, region_end, truth_read, is_kept = tuple(region)

                if not is_kept:
                    continue

                ref_start = region_start
                ref_end = region_end + 1
                # ref_seq should contain region_end_position base
                ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                    ref_start,
                                                                    ref_end)

                read_start = max(0, region_start)
                read_end = region_end
                include_supplementary = False
                all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                       read_start,
                                                       read_end,
                                                       include_supplementary,
                                                       0,
                                                       0)
                total_reads = len(all_reads)

                if total_reads == 0:
                    continue

                total_reads = len(all_reads)
                # sys.stderr.write("INFO: " + log_prefix + " TOTAL " + str(total_reads) + " READS BEFORE DOWNSAMPLING.\n")
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
                # sys.stderr.write("INFO: " + log_prefix + " TOTAL " + str(total_reads) + " READS AFTER DOWNSAMPLING.\n")

                start_time = time.time()

                if realignment_flag:
                    all_reads = self.reads_to_reference_realignment(read_start,
                                                                    read_end,
                                                                    all_reads)
                    # sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " REALIGNMENT OF TOTAL "
                    #                  + str(total_reads) + " READS TOOK: " + str(round(time.time()-start_time, 5))
                    #                  + " secs\n" + TextColor.END)

                summary_generator = PEPPER.SummaryGenerator(ref_seq,
                                                            self.chromosome_name,
                                                            ref_start,
                                                            ref_end)

                summary_generator.generate_train_summary(all_reads,
                                                         region_start,
                                                         region_end,
                                                         truth_read)

                images, labels, positions, chunk_ids = self.chunk_images_train(summary_generator,
                                                                               chunk_size=ImageSizeOptions.SEQ_LENGTH,
                                                                               chunk_overlap=ImageSizeOptions.SEQ_OVERLAP)

                all_images.extend(images)
                all_labels.extend(labels)
                all_positions.extend(positions)
                all_image_chunk_ids.extend(chunk_ids)
        else:
            # HERE REALIGN THE READS TO THE REFERENCE THEN GENERATE THE SUMMARY TO GET A POLISHED HAPLOTYPE
            read_start = max(0, self.region_start_position)
            read_end = self.region_end_position
            include_supplementary = False
            all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                   read_start,
                                                   read_end,
                                                   include_supplementary,
                                                   0,
                                                   0)
            total_reads = len(all_reads)

            if total_reads == 0:
                return [], [], [], []

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

            summary_generator = PEPPER.SummaryGenerator(ref_seq,
                                                        self.chromosome_name,
                                                        self.region_start_position,
                                                        self.region_end_position)

            summary_generator.generate_summary(all_reads,
                                               self.region_start_position,
                                               self.region_end_position)

            images, labels, positions, chunk_ids = \
                self.chunk_images(summary_generator,
                                  chunk_size=ImageSizeOptions.SEQ_LENGTH,
                                  chunk_overlap=ImageSizeOptions.SEQ_OVERLAP)

            all_images.extend(images)
            all_labels.extend(labels)
            all_positions.extend(positions)
            all_image_chunk_ids.extend(chunk_ids)

        # assert(len(all_images) == len(all_labels) == len(all_image_chunk_ids))

        return all_images, all_labels, all_positions, all_image_chunk_ids
