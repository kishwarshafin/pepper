import argparse
import math
import time
import os
import sys
import pickle
import h5py
import numpy as np

from build import HELEN
from modules.python.IntervalTree import IntervalTree
from modules.python.LocalRealignment import LocalAssembler
from modules.python.CandidateFinder import CandidateFinder
from modules.python.TextColor import TextColor
from modules.python.TsvHandler import TsvHandler
from modules.python.FileManager import FileManager
from modules.python.PileupGenerator import PileupGenerator
from modules.python.Options import ImageSizeOptions
from modules.python.CandidateLabler import CandidateLabeler
from modules.python.AlignmentSummarizer import AlignmentSummarizer
"""
This script creates training images from BAM, Reference FASTA and truth VCF file. The process is:
- Find candidates that can be variants
- Label candidates using the VCF
- Create images for each candidate

Input:
- BAM file: Alignment of a genome
- REF file: The reference FASTA file used in the alignment
- VCF file: A truth VCF file
- BED file: A confident bed file. If confident_bed is passed it will only generate train set for those region.

Output:
- TENSOR files: Containing images and their labels.
- CSV file: Containing records of images and their location in the tensor file.
"""

# Global debug helpers
DEBUG_PRINT_CANDIDATES = False
DEBUG_TIME_PROFILE = False
DEBUG_TEST_PARALLEL = False
BED_POSITION_BUFFER = 0


class View:
    """
    Process manager that runs sequence of processes to generate images and their labebls.
    """
    def __init__(self, chromosome_name, bam_file_path, draft_file_path, truth_bam, train_mode, downsample_rate):
        """
        Initialize a manager object
        :param chromosome_name: Name of the chromosome
        :param bam_file_path: Path to the BAM file
        :param draft_file_path: Path to the reference FASTA file
        :param truth_bam: Path to the truth sequence to draft mapping file
        """
        # --- initialize handlers ---
        # create objects to handle different files and query
        self.bam_path = bam_file_path
        self.fasta_path = draft_file_path
        self.truth_bam_path = truth_bam
        self.bam_handler = HELEN.BAM_handler(bam_file_path)
        self.fasta_handler = HELEN.FASTA_handler(draft_file_path)
        self.train_mode = train_mode
        self.downsample_rate = downsample_rate
        self.truth_bam_handler = None
        if self.train_mode:
            self.truth_bam_handler = HELEN.BAM_handler(truth_bam)

        # --- initialize names ---
        # name of the chromosome
        self.chromosome_name = chromosome_name

    @staticmethod
    def build_chromosomal_interval_trees(confident_bed_path):
        """
        Produce a dictionary of intervals trees, with one tree per chromosome
        :param confident_bed_path: Path to confident bed file
        :return: trees_chromosomal
        """
        # create an object for tsv file handling
        tsv_handler_reference = TsvHandler(tsv_file_path=confident_bed_path)
        # create intervals based on chromosome
        intervals_chromosomal_reference = tsv_handler_reference.get_bed_intervals_by_chromosome(universal_offset=-1)

        return intervals_chromosomal_reference

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    @staticmethod
    def a_fully_contains_range_b(range_a, range_b):
        if range_b[0] >= range_a[0] and range_b[1] <= range_a[1]: return True
        return False

    def parse_region(self, start_position, end_position):
        """
        Generate labeled images of a given region of the genome
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :return:
        """
        # 1) First go through all the reads and create a summary
        #   we need: i) Summary of each genomic position
        #           ii) position to index and index to position mapping
        #          iii) easier way to go back and forth
        alignment_summarizer = AlignmentSummarizer(self.bam_handler,
                                                   self.fasta_handler,
                                                   self.chromosome_name,
                                                   start_position,
                                                   end_position)
        images, lables, positions = alignment_summarizer.create_summary(self.truth_bam_handler,
                                                                        self.train_mode)
        return images, lables, positions


def create_output_dir_for_chromosome(output_dir, chr_name):
    """
    Create an internal directory inside the output directory to dump choromosomal summary files
    :param output_dir: Path to output directory
    :param chr_name: chromosome name
    :return: New directory path
    """
    path_to_dir = output_dir + chr_name + "/"
    if not os.path.exists(path_to_dir):
        os.mkdir(path_to_dir)

    summary_path = path_to_dir + "summary" + "/"
    if not os.path.exists(summary_path):
        os.mkdir(summary_path)

    return path_to_dir


def chromosome_level_parallelization(chr_list,
                                     bam_file,
                                     draft_file,
                                     truth_bam,
                                     output_path,
                                     image_path,
                                     total_threads,
                                     thread_id,
                                     train_mode,
                                     downsample_rate,
                                     local_alignment,
                                     log_dir,
                                     max_size=10000):
    """
    This method takes one chromosome name as parameter and chunks that chromosome in max_threads.
    :param chr_list: List of chromosomes to be processed
    :param bam_file: path to BAM file
    :param ref_file: path to reference FASTA file
    :param vcf_file: path to VCF file
    :param max_size: Maximum size of a segment
    :param output_path: path to output directory
    :return:
    """
    # if there's no confident bed provided, then chop the chromosome
    fasta_handler = HELEN.FASTA_handler(draft_file)

    smry = open(output_path + "summary_" + str(thread_id) + ".csv", 'w')

    for chr_name, region in chr_list:
        if not region:
            interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(chr_name) + 1)
        else:
            interval_start, interval_end = tuple(region)

        all_intervals = []
        for pos in range(interval_start, interval_end, max_size):
            all_intervals.append((pos, min(interval_end, pos + max_size - 1)))

        intervals = [r for i, r in enumerate(all_intervals) if i % total_threads == thread_id]

        view = View(chromosome_name=chr_name,
                    bam_file_path=bam_file,
                    draft_file_path=draft_file,
                    truth_bam=truth_bam,
                    train_mode=train_mode,
                    downsample_rate=downsample_rate)

        all_images = []
        all_labels = []
        all_positions = []
        all_indices = []
        global_index = 0
        image_file_name = image_path + chr_name + "_" + str(thread_id) + ".h5py"
        start_time = time.time()

        for count, interval in enumerate(intervals):
            _start, _end = interval
            images, lables, positions = view.parse_region(start_position=_start, end_position=_end)

            # save the images
            for i, image in enumerate(images):
                all_images.append(image)
                all_labels.append(lables[i])
                position_list, index_list = zip(*positions[i])
                all_positions.append(position_list)
                all_indices.append(index_list)

                assert(len(lables[i]) == len(position_list) == len(index_list))
                # write in summary file
                summary_string = image_file_name + "," + str(global_index) + "," + chr_name + "\n"
                smry.write(summary_string)
                global_index += 1

        assert(len(all_images) == len(all_labels) == len(all_positions) == len(all_indices))

        with h5py.File(image_file_name, mode='w') as hdf5_file:
            # the image dataset we save. The index name in h5py is "images".
            img_dset = hdf5_file.create_dataset("image", (len(all_images),) + (ImageSizeOptions.SEQ_LENGTH,
                                                                               ImageSizeOptions.IMAGE_HEIGHT),
                                                np.double, compression='gzip')
            label_dataset = hdf5_file.create_dataset("label", (len(all_labels),) + (ImageSizeOptions.LABEL_LENGTH,),
                                                     np.int16, compression='gzip')
            position_dataset = hdf5_file.create_dataset("position", (len(all_positions),) +
                                                        (ImageSizeOptions.SEQ_LENGTH,), np.int64, compression='gzip')
            index_dataset = hdf5_file.create_dataset("index", (len(all_indices),) + (ImageSizeOptions.SEQ_LENGTH,),
                                                     np.int32, compression='gzip')

            # save the images and labels to the h5py file
            img_dset[...] = all_images
            label_dataset[...] = all_labels
            position_dataset[...] = all_positions
            index_dataset[...] = all_indices

        print("DONE CHROMOSOME: ", chr_name,
              "THREAD ID: ", thread_id,
              "TOTAL TIME ELAPSED: ", int(math.floor(time.time()-start_time)/60), "MINS",
              math.ceil(time.time()-start_time) % 60, "SEC")


def summary_file_to_csv(output_dir_path, chr_list):
    """
    Remove the abundant number of summary files and bind them to one
    :param output_dir_path: Path to the output directory
    :param chr_list: List of chromosomes
    :return:
    """
    for chr_name in chr_list:
        # here we dumped all the bed files
        path_to_dir = output_dir_path + chr_name + "/summary/"

        concatenated_file_name = output_dir_path + chr_name + ".csv"

        filemanager_object = FileManager()
        # get all bed file paths from the directory
        file_paths = filemanager_object.get_file_paths_from_directory(path_to_dir)
        # dump all bed files into one
        filemanager_object.concatenate_files(file_paths, concatenated_file_name)
        # delete all temporary files
        filemanager_object.delete_files(file_paths)
        # remove the directory
        os.rmdir(path_to_dir)


def handle_output_directory(output_dir, thread_id):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    # process the output directory
    if output_dir[-1] != "/":
        output_dir += "/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    internal_directory = "images_" + str(thread_id) + "/"
    image_dir = output_dir + internal_directory

    if not os.path.exists(image_dir):
        os.makedirs(image_dir)

    return output_dir, image_dir


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s:
    :return:
    """
    if s.lower() not in {'false', 'true'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true'


def get_chromosme_list(chromosome_names, ref_file):
    if not chromosome_names:
        fasta_handler = HELEN.FASTA_handler(ref_file)
        chromosome_names = fasta_handler.get_chromosome_names()
        chromosome_names = ','.join(chromosome_names)

    split_names = chromosome_names.strip().split(',')
    split_names = [name.strip() for name in split_names]


    chromosome_name_list = []
    for name in split_names:
        # split on region
        region = None
        if ':' in name:
            name_region = name.strip().split(':')

            if len(name_region) != 2:
                sys.stderr.print(TextColor.RED + "ERROR: --chromosome_name INVALID value.\n" + TextColor.END)
                exit(0)

            name, region = tuple(name_region)
            region = region.strip().split('-')
            region = [int(pos) for pos in region]

            if len(region) != 2 or not region[0] <= region[1]:
                sys.stderr.print(TextColor.RED + "ERROR: --chromosome_name INVALID value.\n" + TextColor.END)
                exit(0)

        range_split = name.split('-')
        if len(range_split) > 1:
            chr_prefix = ''
            for p in name:
                if p.isdigit():
                    break
                else:
                    chr_prefix = chr_prefix + p

            int_ranges = []
            for item in range_split:
                s = ''.join(i for i in item if i.isdigit())
                int_ranges.append(int(s))
            int_ranges = sorted(int_ranges)

            for chr_seq in range(int_ranges[0], int_ranges[-1] + 1):
                chromosome_name_list.append((chr_prefix + str(chr_seq), region))
        else:
            chromosome_name_list.append((name, region))

    return chromosome_name_list


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file containing reads mapped to the draft assembly."
    )
    parser.add_argument(
        "--draft",
        type=str,
        required=True,
        help="FASTA file containing draft assembly."
    )
    parser.add_argument(
        "--truth_bam",
        type=str,
        default=None,
        help="BAM file containing mapping of true reference to the draft assembly"
    )
    parser.add_argument(
        "--chromosome_name",
        type=str,
        help="Desired chromosome number E.g.: 3"
    )
    parser.add_argument(
        "--train_mode",
        type=boolean_string,
        default=False,
        help="If true then a dry test is run."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="candidate_finder_output/",
        help="Path to output directory."
    )
    parser.add_argument(
        "--log_dir",
        type=str,
        default="candidate_finder_output/",
        help="Path to output directory."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )
    parser.add_argument(
        "--thread_id",
        type=int,
        required=False,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--downsample_rate",
        type=float,
        required=False,
        default=1.0,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--local_alignment",
        type=boolean_string,
        default=True,
        help="If true then perform local alignment."
    )
    FLAGS, unparsed = parser.parse_known_args()
    chr_list = get_chromosme_list(FLAGS.chromosome_name, FLAGS.draft)

    if FLAGS.train_mode and (not FLAGS.truth_bam):
        sys.stderr.write(TextColor.RED + "ERROR: TRAIN MODE REQUIRES --vcf AND --bed TO BE SET.\n" + TextColor.END)
        exit(1)
    output_dir, image_dir = handle_output_directory(os.path.abspath(FLAGS.output_dir), FLAGS.thread_id)

    log_dir = FLAGS.log_dir
    if log_dir[-1] != "/":
        log_dir += "/"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    chromosome_level_parallelization(chr_list,
                                     FLAGS.bam,
                                     FLAGS.draft,
                                     FLAGS.truth_bam,
                                     output_dir,
                                     image_dir,
                                     FLAGS.threads,
                                     FLAGS.thread_id,
                                     FLAGS.train_mode,
                                     FLAGS.downsample_rate,
                                     FLAGS.local_alignment,
                                     log_dir)

