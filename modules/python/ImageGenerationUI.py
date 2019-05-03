import math
import time
import os
import sys

from build import HELEN
from modules.python.TextColor import TextColor
from modules.python.DataStore import DataStore
from modules.python.AlignmentSummarizer import AlignmentSummarizer
from modules.python.Options import ImageSizeOptions


class UserInterfaceView:
    """
    Process manager that runs sequence of processes to generate images and their labels.
    """
    def __init__(self, chromosome_name, bam_file_path, draft_file_path, truth_bam, train_mode, downsample_rate):
        """
        Initialize a manager object
        :param chromosome_name: Name of the chromosome
        :param bam_file_path: Path to the BAM file
        :param draft_file_path: Path to the reference FASTA file
        :param truth_bam: Path to the truth sequence to reference mapping file
        """
        # --- initialize handlers ---
        # create objects to handle different files and query
        self.bam_path = bam_file_path
        self.fasta_path = draft_file_path
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

    def parse_region(self, start_position, end_position):
        """
        Generate labeled images of a given region of the genome
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :return:
        """
        alignment_summarizer = AlignmentSummarizer(self.bam_handler,
                                                   self.fasta_handler,
                                                   self.chromosome_name,
                                                   start_position,
                                                   end_position)

        images, lables, positions, image_chunk_ids = alignment_summarizer.create_summary(self.truth_bam_handler,
                                                                                         self.train_mode,
                                                                                         realignment_flag=False)

        return images, lables, positions, image_chunk_ids


class UserInterfaceSupport:
    """
    THIS CLASS HOLDS STATIC METHODS THAT HELP AS USER INTERFACE FOR IMAGE GENERATION
    """
    @staticmethod
    def chunks(intervals, max_chunk_size):
        """Yield successive n-sized chunks from l."""
        chunks = []
        for i in range(0, len(intervals), max_chunk_size):
            chunks.append(intervals[i:i + max_chunk_size])
        return chunks

    @staticmethod
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

    @staticmethod
    def get_chromosome_list(chromosome_names, ref_file):
        """
        PARSES THROUGH THE CHROMOSOME PARAMETER TO FIND OUT WHICH REGIONS TO PROCESS
        :param chromosome_names: NAME OF CHROMOSOME
        :param ref_file: PATH TO THE REFERENCE FILE
        :return: LIST OF CHROMOSOME IN REGION SPECIFIC FORMAT
        """
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

    @staticmethod
    def single_worker(args, _start, _end):
        chr_name, bam_file, draft_file, truth_bam, train_mode, downsample_rate = args

        view = UserInterfaceView(chromosome_name=chr_name,
                                 bam_file_path=bam_file,
                                 draft_file_path=draft_file,
                                 truth_bam=truth_bam,
                                 train_mode=train_mode,
                                 downsample_rate=downsample_rate)

        images, lables, positions, image_chunk_ids = view.parse_region(_start, _end)

        region = (chr_name, _start, _end)

        return images, lables, positions, image_chunk_ids, region

    @staticmethod
    def chromosome_level_parallelization(chr_list,
                                         bam_file,
                                         draft_file,
                                         truth_bam,
                                         output_path,
                                         total_threads,
                                         thread_id,
                                         train_mode,
                                         downsample_rate,
                                         max_size=10000):
        start_time = time.time()
        fasta_handler = HELEN.FASTA_handler(draft_file)

        file_name = output_path + "helen_images_" + str(thread_id) + ".hdf"
        thread_prefix = "{:02d}".format(thread_id) + "/" + "{:02d}".format(total_threads) + ":"

        with DataStore(file_name, 'w') as output_hdf_file:
            for chr_name, region in chr_list:
                # sys.stderr.write(TextColor.CYAN + "INFO: " + thread_prefix + " GENERATING IMAGE FROM CONTIG "
                #                  + str(chr_name) + "\n" + TextColor.END)
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
                    pos_start = max(interval_start, pos - ImageSizeOptions.MIN_IMAGE_OVERLAP)
                    pos_end = min(interval_end, pos + max_size + ImageSizeOptions.MIN_IMAGE_OVERLAP)
                    all_intervals.append((pos_start, pos_end))

                args = (chr_name, bam_file, draft_file, truth_bam, train_mode, downsample_rate)

                intervals = [r for i, r in enumerate(all_intervals) if i % total_threads == thread_id]

                for _start, _end in intervals:
                    images, labels, positions, chunk_ids, region = UserInterfaceSupport.single_worker(args, _start, _end)

                    log_prefix = "[" + str(region[0]) + ":" + str(region[1]) + "-" + str(region[2]) + "]"
                    for i, image in enumerate(images):
                        label = labels[i]
                        position, index = zip(*positions[i])
                        chunk_id = chunk_ids[i]

                        summary_name = str(region[0]) + "_" + str(region[1]) + "_" + str(region[2]) + "_" + str(chunk_id)

                        output_hdf_file.write_summary(region, image, label, position, index, chunk_id, summary_name)

                    # sys.stderr.write(TextColor.GREEN + "INFO: " + thread_prefix + " " + log_prefix + " TOTAL "
                    #                  + str(len(images)) + " IMAGES SAVED\n" + TextColor.END)

                # sys.stderr.write(TextColor.BLUE + "INFO: " + thread_prefix + " COMPLETED PROCESSING CHROMOSOME: " +
                #                  chr_name + " TOTAL TIME ELAPSED: " + str(int(math.floor(time.time()-start_time)/60))
                #                  + " MINS " + str(math.ceil(time.time()-start_time) % 60) + " SEC\n" + TextColor.END)
