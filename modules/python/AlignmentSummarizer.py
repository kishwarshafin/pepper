from build import HELEN
import numpy as np
import math


class AlignmentSummarizer:
    def __init__(self, bam_handler, fasta_handler, chromosome_name, region_start, region_end):
        self.bam_handler = bam_handler
        self.fasta_handler = fasta_handler
        self.chromosome_name = chromosome_name
        self.region_start_position = region_start
        self.region_end_position = region_end

    def chunk_images(self, summary, train_mode, chunk_size, chunk_overlap):

        chunk_start = 0
        chunk_end = min(len(summary.genomic_pos), chunk_size)
        images = []
        labels = []
        positions = []

        while True:
            image_chunk = summary.image[chunk_start:chunk_end]
            pos_chunk = summary.genomic_pos[chunk_start:chunk_end]
            if train_mode:
                label_chunk = summary.labels[chunk_start:chunk_end]
            else:
                label_chunk = [0] * (chunk_end-chunk_start)

            assert(len(image_chunk) == len(pos_chunk) == len(label_chunk))
            # print(len(image_chunk), len(pos_chunk), len(label_chunk))

            padding_required = chunk_size - len(image_chunk)
            if padding_required > 0:
                label_chunk = label_chunk + [0] * padding_required
                pos_chunk = pos_chunk + [(0, 0)] * padding_required
                image_chunk = image_chunk + [[0.0] * 20] * padding_required

            assert(len(image_chunk) == len(pos_chunk) == len(label_chunk))

            if chunk_end == len(summary.genomic_pos):
                break

            chunk_start = chunk_end - chunk_overlap
            chunk_end = min(len(summary.genomic_pos), chunk_start + chunk_size)
            images.append(image_chunk)
            labels.append(label_chunk)
            positions.append(pos_chunk)

        return images, labels, positions

    def create_summary(self, truth_bam_handler, train_mode):
        truth_reads = []
        if train_mode:
            truth_reads = truth_bam_handler.get_reads(self.chromosome_name,
                                                      self.region_start_position,
                                                      self.region_end_position,
                                                      0,
                                                      0)
            if len(truth_reads) > 1:
                print("MORE THAN ONE TRUTH READ IN REGION: ", self.chromosome_name,
                      self.region_start_position,
                      self.region_end_position)

        # get the reads from the bam file
        all_reads = self.bam_handler.get_reads(self.chromosome_name,
                                               self.region_start_position,
                                               self.region_end_position,
                                               0,
                                               0)

        # ref_seq should contain region_end_position base
        ref_seq = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                            self.region_start_position,
                                                            self.region_end_position + 1)

        summary_generator = HELEN.SummaryGenerator(ref_seq,
                                                   self.chromosome_name,
                                                   self.region_start_position,
                                                   self.region_end_position)
        summary_generator.generate_summary(all_reads,
                                           self.region_start_position,
                                           self.region_end_position,
                                           truth_reads,
                                           train_mode)

        images, labels, positions = self.chunk_images(summary_generator, train_mode, chunk_size=1000, chunk_overlap=200)

        return images, labels, positions
