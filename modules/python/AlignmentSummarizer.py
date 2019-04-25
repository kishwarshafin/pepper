from build import HELEN
import itertools
from operator import itemgetter
import sys
from modules.python.TextColor import TextColor
from modules.python.Options import ImageSizeOptions


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
        chunk_end = min(len(summary.genomic_pos), chunk_size)
        images = []
        labels = []
        positions = []

        while True:
            image_chunk = summary.image[chunk_start:chunk_end]
            pos_chunk = summary.genomic_pos[chunk_start:chunk_end]
            label_chunk = [0] * (chunk_end - chunk_start)

            assert (len(image_chunk) == len(pos_chunk) == len(label_chunk))
            # print(len(image_chunk), len(pos_chunk), len(label_chunk))

            padding_required = chunk_size - len(image_chunk)
            if padding_required > 0:
                label_chunk = label_chunk + [0] * padding_required
                pos_chunk = pos_chunk + [(-1, -1)] * padding_required
                image_chunk = image_chunk + [[0.0] * ImageSizeOptions.IMAGE_HEIGHT] * padding_required

            assert (len(image_chunk) == len(pos_chunk) == len(label_chunk) == ImageSizeOptions.SEQ_LENGTH)

            images.append(image_chunk)
            labels.append(label_chunk)
            positions.append(pos_chunk)

            if chunk_end == len(summary.genomic_pos):
                break

            chunk_start = chunk_end - chunk_overlap
            chunk_end = min(len(summary.genomic_pos), chunk_start + chunk_size)

        return images, labels, positions

    @staticmethod
    def chunk_images_train(summary, chunk_size, chunk_overlap):
        chunk_start = 0
        chunk_end = min(len(summary.genomic_pos), chunk_size)
        images = []
        labels = []
        positions = []

        while True:
            if chunk_end - chunk_start != chunk_size:
                padding_required = chunk_size - (chunk_end - chunk_start)
                chunk_start -= padding_required
                if chunk_start < 0:
                    break

            image_chunk = summary.image[chunk_start:chunk_end]
            pos_chunk = summary.genomic_pos[chunk_start:chunk_end]
            label_chunk = summary.labels[chunk_start:chunk_end]

            assert (len(image_chunk) == len(pos_chunk) == len(label_chunk) == chunk_size)

            images.append(image_chunk)
            labels.append(label_chunk)
            positions.append(pos_chunk)

            if chunk_end == len(summary.genomic_pos):
                break

            chunk_start = chunk_end - chunk_overlap
            chunk_end = min(len(summary.genomic_pos), chunk_start + chunk_size)

        return images, labels, positions

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    def get_overlap_between_ranges(self, range_a, range_b):
        if range_a[1] > range_b[0]:
            return range_b[0], range_a[1]
        else:
            return None

    def remove_conflicting_regions(self, regions, min_length=1000, length_ratio=2.0, overlap_fraction=0.5):
        # reused from medaka's filter_alignments method.
        for reg_a, reg_b in itertools.combinations(regions, 2):
            el1, el2 = sorted((reg_a, reg_b), key=itemgetter(0))
            overlap = self.get_overlap_between_ranges(el1, el2)

            if overlap is None:
                continue
            ovlp_start, ovlp_end = overlap
            s, l = sorted((reg_a, reg_b), key=lambda element: (element[1] - element[0]))

            length_ratio_ij = (l[1] - l[0]) / (s[1] - s[0])
            overlap_fraction_ij = (ovlp_end - ovlp_start) / (s[1] - s[0])
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

    def create_summary(self, truth_bam_handler, train_mode):
        all_images = []
        all_labels = []
        all_positions = []
        if train_mode:
            truth_reads = truth_bam_handler.get_reads(self.chromosome_name,
                                                      self.region_start_position,
                                                      self.region_end_position,
                                                      0,
                                                      0)
            regions = []
            for read in truth_reads:
                regions.append([read.pos, read.pos_end, read, True])

            regions = self.remove_conflicting_regions(regions, min_length=ImageSizeOptions.SEQ_LENGTH)

            for region in regions:
                region_start, region_end, truth_read, is_kept = tuple(region)

                if not is_kept:
                    continue

                # get the reads from the bam file
                all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                       region_start,
                                                       region_end,
                                                       0,
                                                       0)

                # ref_seq should contain region_end_position base
                ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                    region_start,
                                                                    region_end + 1)

                summary_generator = HELEN.SummaryGenerator(ref_seq,
                                                           self.chromosome_name,
                                                           region_start,
                                                           region_end)

                summary_generator.generate_train_summary(all_reads,
                                                         region_start,
                                                         region_end,
                                                         truth_read,
                                                         train_mode)

                images, labels, positions = self.chunk_images_train(summary_generator,
                                                                    chunk_size=ImageSizeOptions.SEQ_LENGTH,
                                                                    chunk_overlap=ImageSizeOptions.SEQ_OVERLAP)

                all_images.extend(images)
                all_labels.extend(labels)
                all_positions.extend(positions)
        else:
            # get the reads from the bam file
            all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                                   self.region_start_position,
                                                   self.region_end_position,
                                                   0,
                                                   0)

            reads_un, reads_hp1, reads_hp2 = all_reads
            total_reads = len(reads_un) + len(reads_hp1) + len(reads_hp2)

            sys.stderr.write(TextColor.GREEN + "INFO: TOTAL " + str(total_reads) + " READS FOUND IN " +
                             str(self.chromosome_name) + ":" + str(self.region_start_position) + "-"
                             + str(self.region_end_position) + " WITH: " + str(len(reads_un)) + " UNTAGGED, "
                             + str(len(reads_hp1)) + " HAPLOTYPE_1, " + str(len(reads_hp2)) + " HAPLOTYPE_2 READS.\n"
                             + TextColor.END)

            # HERE REALIGN THE READS TO THE REFERENCE THEN GENERATE THE SUMMARY TO GET A POLISHED HAPLOTYPE
            # ref_seq should contain region_end_position base
            ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                self.region_start_position,
                                                                self.region_end_position + 1)

            summary_generator = HELEN.SummaryGenerator(ref_seq,
                                                       self.chromosome_name,
                                                       self.region_start_position,
                                                       self.region_end_position)

            summary_generator.generate_summary(reads_un + reads_hp1+ reads_hp2,
                                               self.region_start_position,
                                               self.region_end_position)
            exit()

            images, labels, positions = self.chunk_images(summary_generator,
                                                          chunk_size=ImageSizeOptions.SEQ_LENGTH,
                                                          chunk_overlap=ImageSizeOptions.SEQ_OVERLAP)

            all_images.extend(images)
            all_labels.extend(labels)
            all_positions.extend(positions)

        return all_images, all_labels, all_positions
