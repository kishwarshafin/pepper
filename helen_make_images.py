import argparse
import math
import time
import os
import sys

from build import HELEN
from modules.python.IntervalTree import IntervalTree
from modules.python.TextColor import TextColor
from modules.python.TsvHandler import TsvHandler
from modules.python.Options import ImageSizeOptions
from modules.python.AlignmentSummarizer import AlignmentSummarizer
from modules.python.datastore import DataStore
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
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
    def __init__(self, chromosome_name, bam_file_path, draft_file_path, truth_bam_h1, truth_bam_h2,
                 train_mode, downsample_rate):
        """
        Initialize a manager object
        :param chromosome_name: Name of the chromosome
        :param bam_file_path: Path to the BAM file
        :param draft_file_path: Path to the reference FASTA file
        :param truth_bam_h1: Path to the truth sequence to reference mapping file
        :param truth_bam_h2: Path to the truth sequence to reference mapping file
        """
        # --- initialize handlers ---
        # create objects to handle different files and query
        self.bam_path = bam_file_path
        self.fasta_path = draft_file_path
        self.bam_handler = HELEN.BAM_handler(bam_file_path)
        self.fasta_handler = HELEN.FASTA_handler(draft_file_path)
        self.train_mode = train_mode
        self.downsample_rate = downsample_rate
        self.truth_bam_handler_h1 = None
        self.truth_bam_handler_h2 = None
        if self.train_mode:
            self.truth_bam_handler_h1 = HELEN.BAM_handler(truth_bam_h1)
            self.truth_bam_handler_h2 = HELEN.BAM_handler(truth_bam_h2)

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

        images, lables, positions, image_hp_tags, image_chunk_ids = \
            alignment_summarizer.create_summary(self.truth_bam_handler_h1,
                                                self.truth_bam_handler_h2,
                                                self.train_mode,
                                                realignment_flag=False)

        return images, lables, positions, image_hp_tags, image_chunk_ids


def single_worker(args, _start, _end, confident_bed_regions):
    chr_name, bam_file, draft_file, truth_bam_h1, truth_bam_h2, train_mode, downsample_rate = args
    if train_mode:
        # if its in train mode then split it into intervals of confident regions
        intervals = get_confident_intervals_of_a_region(confident_bed_regions, _start, _end)
    else:
        # otherwise simply chunk it using the same logic as before but the chunk size is 10^4 this time
        max_size = 10000
        intervals = []
        for pos in range(_start, _end, max_size):
            intervals.append((pos, min(_end, pos + max_size)))

    all_images = []
    all_labels = []
    all_positions = []
    all_hp_tags = []
    all_chunk_ids = []
    for region_start, region_end in intervals:
        view = View(chromosome_name=chr_name,
                    bam_file_path=bam_file,
                    draft_file_path=draft_file,
                    truth_bam_h1=truth_bam_h1,
                    truth_bam_h2=truth_bam_h2,
                    train_mode=train_mode,
                    downsample_rate=downsample_rate)
        images, lables, positions, image_hp_tags, image_chunk_ids = view.parse_region(region_start, region_end)

        all_images.extend(images)
        all_labels.extend(lables)
        all_positions.extend(positions)
        all_hp_tags.extend(image_hp_tags)
        all_chunk_ids.extend(image_chunk_ids)

    region = (chr_name, _start, _end)

    return all_images, all_labels, all_positions, all_hp_tags, region, all_chunk_ids


def get_confident_intervals_of_a_region(confident_bed_regions, start_position, end_position,
                                        min_size=ImageSizeOptions.MIN_SEQUENCE_LENGTH, max_size=10000):

    confident_intervals_in_region = confident_bed_regions.find(start_position, end_position)

    all_intervals = []
    for interval_start, interval_end in confident_intervals_in_region:
        if interval_end < start_position:
            continue
        if interval_start > end_position:
            continue

        interval_start = max(interval_start, start_position)
        interval_end = min(interval_end, end_position)

        if interval_end - interval_start > max_size:
            for pos in range(interval_start, interval_end, max_size):
                all_intervals.append((pos, min(interval_end, pos + max_size)))
        elif interval_end - interval_start >= min_size:
            all_intervals.append((interval_start, interval_end))
    return all_intervals


def chunks(intervals, max_chunk_size):
    """Yield successive n-sized chunks from l."""
    chunks = []
    for i in range(0, len(intervals), max_chunk_size):
        chunks.append(intervals[i:i + max_chunk_size])
    return chunks


def chromosome_level_parallelization(chr_list,
                                     bam_file,
                                     draft_file,
                                     truth_bam_h1,
                                     truth_bam_h2,
                                     confident_bed_regions,
                                     output_path,
                                     total_threads,
                                     train_mode,
                                     downsample_rate,
                                     max_size=100000):
    start_time = time.time()
    fasta_handler = HELEN.FASTA_handler(draft_file)
    chr_counter = 1

    timestr = time.strftime("%m%d%Y_%H%M%S")
    file_name = output_path + "helen_images_" + timestr + ".hdf"

    with DataStore(file_name, 'w') as output_hdf_file:
        for chr_name, region in chr_list:
            sys.stderr.write(TextColor.GREEN + "INFO: " + str(chr_counter) + "/" + str(len(chr_list))
                             + " GENERATING IMAGE FROM CONTIG " + str(chr_name) + "\n" + TextColor.END)
            if not region:
                interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(chr_name) - 1)
                max_end = interval_end
            else:
                interval_start, interval_end = tuple(region)
                interval_start = max(0, interval_start)
                interval_end = min(interval_end, fasta_handler.get_chromosome_sequence_length(chr_name) - 1)
                max_end = interval_end

            # this is the interval size each of the process is going to get which is 10^6
            # I will split this into 10^4 size inside the worker process
            all_intervals = []
            for pos in range(interval_start, interval_end, max_size):
                all_intervals.append((pos, min(interval_end, pos + max_size)))

            confident_tree = None
            if train_mode:
                confident_tree = IntervalTree(confident_bed_regions[chr_name])

            args = (chr_name, bam_file, draft_file, truth_bam_h1, truth_bam_h2, train_mode, downsample_rate)

            with ProcessPoolExecutor(max_workers=total_threads) as executor:
                futures = [executor.submit(single_worker, args, _start, _end, confident_tree) for _start, _end in
                           all_intervals]

                for fut in as_completed(futures):
                    if fut.exception() is None:
                        images, labels, positions, image_hp_tag, region, chunk_ids = fut.result()
                        log_prefix = "[" + str(region[0]) + ":" + str(region[2]) + "/" + str(max_end) + "]"
                        sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " TOTAL " + str(len(images))
                                         + " IMAGES GENERATED\n" + TextColor.END)

                        for i, image in enumerate(images):
                            label = labels[i]
                            position, index = zip(*positions[i])
                            hp_tag = image_hp_tag[i]
                            chunk_id = chunk_ids[i]

                            summary_name = str(region[0]) + "_" + str(region[1]) + "_" + str(region[2]) + "_" + \
                                           str(hp_tag) + "_" + str(i)

                            output_hdf_file.write_summary(region, image, label, position, index, hp_tag, chunk_id,
                                                          summary_name)
                        sys.stderr.write(TextColor.GREEN + "INFO: " + log_prefix + " IMAGES SAVED\n" + TextColor.END)
                        del images, labels, positions, image_hp_tag, region, chunk_ids
                    else:
                        sys.stderr.write(TextColor.RED + "EXCEPTION: " + str(fut.exception()) + "\n" + TextColor.END)
                    fut._result = None

            print("DONE CHROMOSOME: ", chr_name,
                  "TOTAL TIME ELAPSED: ", int(math.floor(time.time()-start_time)/60), "MINS",
                  math.ceil(time.time()-start_time) % 60, "SEC")


def handle_output_directory(output_directory):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_directory: Output directory path
    :return:
    """
    # process the output directory
    if output_directory[-1] != "/":
        output_directory += "/"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    return output_directory


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s:
    :return:
    """
    if s.lower() not in {'false', 'true'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true'


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


def get_chromosome_list(chromosome_names, ref_file):
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
        help="BAM file containing haplotyped (With HP tag) reads mapped to the reference."
    )
    parser.add_argument(
        "--reference",
        type=str,
        required=True,
        help="FASTA file containing the reference sequence."
    )
    parser.add_argument(
        "--truth_bam_h1",
        type=str,
        default=None,
        help="BAM file containing mapping of true haplotype 1 to the reference"
    )
    parser.add_argument(
        "--truth_bam_h2",
        type=str,
        default=None,
        help="BAM file containing mapping of true haplotype 2 to the reference"
    )
    parser.add_argument(
        "--bed",
        type=str,
        default=None,
        help="BED file containing confident regions."
    )
    parser.add_argument(
        "--chromosome_name",
        type=str,
        help="Desired chromosome number [chr_name:start-end] E.g.: chr3:1000-2000 or simply chr3"
    )
    parser.add_argument(
        "--train_mode",
        type=boolean_string,
        default=False,
        help="Generate training images"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="candidate_finder_output/",
        help="Path to output directory."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=5,
        help="Number of maximum threads to use."
    )
    parser.add_argument(
        "--downsample_rate",
        type=float,
        required=False,
        default=1.0,
        help="Downsample reads by this margin."
    )
    FLAGS, unparsed = parser.parse_known_args()
    chr_list = get_chromosome_list(FLAGS.chromosome_name, FLAGS.reference)

    confident_regions = []
    if FLAGS.train_mode:
        confident_regions = build_chromosomal_interval_trees(FLAGS.bed)
        if not FLAGS.truth_bam_h1 or not FLAGS.truth_bam_h2 or not FLAGS.bed:
            raise Exception(TextColor.RED + "ERROR: TRAIN MODE REQUIRES --truth_bam_h1, --truth_bam_h2 "
                                            "AND --bed TO BE SET.\n" + TextColor.END)

    output_dir = handle_output_directory(os.path.abspath(FLAGS.output_dir))

    chromosome_level_parallelization(chr_list,
                                     FLAGS.bam,
                                     FLAGS.reference,
                                     FLAGS.truth_bam_h1,
                                     FLAGS.truth_bam_h2,
                                     confident_regions,
                                     output_dir,
                                     FLAGS.threads,
                                     FLAGS.train_mode,
                                     FLAGS.downsample_rate)

