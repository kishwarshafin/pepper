import time
import os
import sys
import re
import concurrent.futures
from collections import defaultdict
from datetime import datetime
from pepper_snp.build import PEPPER_SNP
from pepper_snp.modules.python.ExcludeContigs import EXCLUDED_HUMAN_CONTIGS
from pepper_snp.modules.python.DataStore import DataStore
from pepper_snp.modules.python.AlignmentSummarizer import AlignmentSummarizer
from pepper_snp.modules.python.Options import ImageSizeOptions
from pepper_snp.modules.python.ExcludeContigs import EXCLUDED_HUMAN_CONTIGS


class UserInterfaceView:
    """
    Process manager that runs sequence of processes to generate images and their labels.
    """
    def __init__(self, chromosome_name, bam_file_path, draft_file_path, truth_bam_h1, truth_bam_h2, realignment_flag,
                 train_mode):
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
        self.bam_handler = PEPPER_SNP.BAM_handler(bam_file_path)
        self.fasta_handler = PEPPER_SNP.FASTA_handler(draft_file_path)
        self.train_mode = train_mode
        self.realignment_flag = realignment_flag
        self.downsample_rate = 1.0
        self.truth_bam_handler_h1 = None
        self.truth_bam_handler_h2 = None

        if self.train_mode:
            self.truth_bam_handler_h1 = PEPPER_SNP.BAM_handler(truth_bam_h1)
            self.truth_bam_handler_h2 = PEPPER_SNP.BAM_handler(truth_bam_h2)

        # --- initialize names ---
        # name of the chromosome
        self.chromosome_name = chromosome_name

    def parse_region(self, start_position, end_position, downsample_rate, bed_list):
        """
        Generate labeled images of a given region of the genome
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :param downsample_rate: Read downsample rate
        :param bed_list: List from bed file
        :return:
        """
        alignment_summarizer = AlignmentSummarizer(self.bam_handler,
                                                   self.fasta_handler,
                                                   self.chromosome_name,
                                                   start_position,
                                                   end_position)

        images, lables, positions, image_chunk_ids, ref_seqs = \
            alignment_summarizer.create_summary(self.truth_bam_handler_h1,
                                                self.truth_bam_handler_h2,
                                                self.train_mode,
                                                downsample_rate=downsample_rate,
                                                bed_list=bed_list,
                                                realignment_flag=self.realignment_flag)

        return images, lables, positions, image_chunk_ids, ref_seqs


class UserInterfaceSupport:
    """
    THIS CLASS HOLDS STATIC METHODS THAT HELP AS USER INTERFACE FOR IMAGE GENERATION
    """
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
    def natural_key(string_):
        """See http://www.codinghorror.com/blog/archives/001018.html"""
        return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

    @staticmethod
    def get_chromosome_list(chromosome_names, ref_file, bam_file, region_bed):
        """
        PARSES THROUGH THE CHROMOSOME PARAMETER TO FIND OUT WHICH REGIONS TO PROCESS
        :param chromosome_names: NAME OF CHROMOSOME
        :param ref_file: PATH TO BAM FILE
        :param bam_file: PATH TO THE REFERENCE FILE
        :param region_bed: PATH TO A BED FILE
        :return: LIST OF CHROMOSOME IN REGION SPECIFIC FORMAT
        """
        chromosome_name_list = []
        region_bed_list = None

        if not chromosome_names:
            fasta_handler = PEPPER_SNP.FASTA_handler(ref_file)
            bam_handler = PEPPER_SNP.BAM_handler(bam_file)
            bam_contigs = bam_handler.get_chromosome_sequence_names()
            fasta_contigs = fasta_handler.get_chromosome_names()
            common_contigs = list(set(fasta_contigs) & set(bam_contigs))
            common_contigs = list(set(common_contigs) - set(EXCLUDED_HUMAN_CONTIGS))

            if len(common_contigs) == 0:
                sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] "
                                 + "ERROR: NO COMMON CONTIGS FOUND BETWEEN THE BAM FILE AND THE FASTA FILE.")
                sys.stderr.flush()
                exit(1)

            common_contigs = sorted(common_contigs, key=UserInterfaceSupport.natural_key)
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: COMMON CONTIGS FOUND: " + str(common_contigs) + "\n")
            sys.stderr.flush()

            for contig_name in common_contigs:
                chromosome_name_list.append((contig_name, None))
        else:
            split_names = chromosome_names.strip().split(',')
            split_names = [name.strip() for name in split_names]
            chromosome_name_list = []
            for name in split_names:
                # split on region
                region = None
                if ':' in name:
                    name_region = name.strip().split(':')

                    if len(name_region) != 2:
                        sys.stderr.write("ERROR: --region INVALID value.\n")
                        exit(0)

                    name, region = tuple(name_region)
                    region = region.strip().split('-')
                    region = [int(pos) for pos in region]

                    if len(region) != 2 or not region[0] <= region[1]:
                        sys.stderr.write("ERROR: --region INVALID value.\n")
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

        if region_bed:
            region_bed_list = defaultdict()
            with open(region_bed) as fp:
                line = fp.readline()
                cnt = 1
                while line:
                    line_to_list = line.rstrip().split('\t')
                    chr_name, start_pos, end_pos = line_to_list[0], int(line_to_list[1]), int(line_to_list[2])
                    region = sorted([start_pos, end_pos])
                    if chr_name not in region_bed_list.keys():
                        region_bed_list[chr_name] = []
                    region_bed_list[chr_name].append(region)
                    line = fp.readline()
                cnt += 1

        return chromosome_name_list, region_bed_list

    @staticmethod
    def single_worker(args, _start, _end):
        chr_name, bam_file, draft_file, truth_bam_h1, truth_bam_h2, realignment_flag, train_mode, downsample_rate, bed_list = args

        view = UserInterfaceView(chromosome_name=chr_name,
                                 bam_file_path=bam_file,
                                 draft_file_path=draft_file,
                                 truth_bam_h1=truth_bam_h1,
                                 truth_bam_h2=truth_bam_h2,
                                 train_mode=train_mode,
                                 realignment_flag=realignment_flag)

        images, labels, positions, image_chunk_ids, ref_seqs = view.parse_region(_start, _end, downsample_rate, bed_list)
        region = (chr_name, _start, _end)

        return images, labels, positions, image_chunk_ids, ref_seqs, region

    @staticmethod
    def image_generator(args, all_intervals, total_threads, thread_id):
        thread_prefix = "[THREAD " + "{:02d}".format(thread_id) + "]"

        output_path, bam_file, draft_file, truth_bam_h1, truth_bam_h2, realignment_flag, train_mode, downsample_rate, bed_list = args
        timestr = time.strftime("%m%d%Y_%H%M%S")
        file_name = output_path + "pepper_snp_images_thread_" + str(thread_id) + "_" + str(timestr) + ".hdf"

        intervals = [r for i, r in enumerate(all_intervals) if i % total_threads == thread_id]

        # initial notification
        if thread_id == 0:
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] "
                             + "STARTING THREAD: " + str(thread_id)
                             + " FOR " + str(len(intervals)) + " INTERVALS\n")
            sys.stderr.flush()

        start_time = time.time()
        with DataStore(file_name, 'w') as output_hdf_file:
            for counter, interval in enumerate(intervals):
                chr_name, _start, _end = interval
                img_args = (chr_name, bam_file, draft_file, truth_bam_h1, truth_bam_h2, realignment_flag, train_mode, downsample_rate, bed_list)
                images, labels, positions, chunk_ids, ref_seqs, region = \
                    UserInterfaceSupport.single_worker(img_args, _start, _end)

                for i, image in enumerate(images):
                    label = labels[i]
                    position, index = zip(*positions[i])
                    chunk_id = chunk_ids[i]
                    ref_seq = ref_seqs[i]

                    summary_name = str(region[0]) + "_" + str(region[1]) + "_" + str(region[2]) + "_" + str(chunk_id)

                    output_hdf_file.write_summary(region, image, label, position, index, chunk_id, ref_seq, summary_name)

                if counter > 0 and counter % 10 == 0:
                    percent_complete = int((100 * counter) / len(intervals))
                    time_now = time.time()
                    mins = int((time_now - start_time) / 60)
                    secs = int((time_now - start_time)) % 60
                    if thread_id == 0:
                        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "]"
                                         + " INFO: " + " " + str(counter) + "/" + str(len(intervals))
                                         + " COMPLETE (" + str(percent_complete) + "%)"
                                         + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                        sys.stderr.flush()

        return thread_id

    @staticmethod
    def chromosome_level_parallelization(chr_list,
                                         bam_file,
                                         draft_file,
                                         truth_bam_h1,
                                         truth_bam_h2,
                                         output_path,
                                         total_threads,
                                         realignment_flag,
                                         train_mode,
                                         downsample_rate,
                                         bed_list):

        if train_mode:
            max_size = 100000
        else:
            max_size = 100000

        start_time = time.time()
        fasta_handler = PEPPER_SNP.FASTA_handler(draft_file)

        all_intervals = []
        # first calculate all the intervals that we need to process
        for chr_name, region in chr_list:
            # contig update message
            if not region:
                interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(chr_name) - 1)
            else:
                interval_start, interval_end = tuple(region)
                interval_start = max(0, interval_start)
                interval_end = min(interval_end, fasta_handler.get_chromosome_sequence_length(chr_name) - 1)

            # this is the interval size each of the process is going to get which is 10^6
            # I will split this into 10^4 size inside the worker process
            for pos in range(interval_start, interval_end, max_size):
                pos_start = max(interval_start, pos - ImageSizeOptions.MIN_IMAGE_OVERLAP)
                pos_end = min(interval_end, pos + max_size + ImageSizeOptions.MIN_IMAGE_OVERLAP)
                all_intervals.append((chr_name, pos_start, pos_end))

        # all intervals calculated now
        # contig update message
        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] "
                         + "INFO: TOTAL CONTIGS: " + str(len(chr_list))
                         + " TOTAL INTERVALS: " + str(len(all_intervals)) + "\n")
        sys.stderr.flush()

        args = (output_path, bam_file, draft_file, truth_bam_h1, truth_bam_h2, realignment_flag, train_mode, downsample_rate, bed_list)
        with concurrent.futures.ProcessPoolExecutor(max_workers=total_threads) as executor:
            futures = [executor.submit(UserInterfaceSupport.image_generator, args, all_intervals, total_threads, thread_id)
                       for thread_id in range(0, total_threads)]

            for fut in concurrent.futures.as_completed(futures):
                if fut.exception() is None:
                    # get the results
                    thread_id = fut.result()
                    if thread_id == 0:
                        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] "
                                         + "THREAD " + str(thread_id) + " FINISHED SUCCESSFULLY.\n")
                else:
                    sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
                fut._result = None  # python issue 27144

        end_time = time.time()
        mins = int((end_time - start_time) / 60)
        secs = int((end_time - start_time)) % 60
        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] FINISHED IMAGE GENERATION\n")
        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] TOTAL ELAPSED TIME FOR IMAGE GENERATION: " + str(mins) + " Min " + str(secs) + " Sec\n")
