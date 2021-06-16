import time
import os
import re
import sys
import pickle
import concurrent.futures
from datetime import datetime
from collections import defaultdict
from pepper_variant.build import PEPPER_VARIANT
from pepper_variant.modules.python.ExcludeContigs import EXCLUDED_HUMAN_CONTIGS
from pepper_variant.modules.python.DataStore import DataStore
from pepper_variant.modules.python.AlignmentSummarizer import AlignmentSummarizer
from pepper_variant.modules.python.AlignmentSummarizerHP import AlignmentSummarizerHP
from pepper_variant.modules.python.Options import ImageSizeOptions


class ImageGenerator:
    """
    Process manager that runs sequence of processes to generate images and their labels.
    """
    def __init__(self, chromosome_name, bam_file_path, draft_file_path, use_hp_info, truth_vcf, train_mode, random_draw_probability):
        """
        Initialize a manager object
        :param chromosome_name: Name of the chromosome
        :param bam_file_path: Path to the BAM file
        :param draft_file_path: Path to the reference FASTA file
        :param truth_bam_h1: Path to the truth sequence of hp1 to reference mapping file
        :param truth_bam_h2: Path to the truth sequence of hp2 to reference mapping file
        """
        # --- initialize handlers ---
        # create objects to handle different files and query
        self.bam_path = bam_file_path
        self.fasta_path = draft_file_path
        self.bam_handler = PEPPER_VARIANT.BAM_handler(bam_file_path)
        self.fasta_handler = PEPPER_VARIANT.FASTA_handler(draft_file_path)
        self.train_mode = train_mode
        self.downsample_rate = 1.0
        self.use_hp_info = use_hp_info
        self.truth_vcf = None
        self.random_draw_probability = random_draw_probability

        if self.train_mode:
            self.truth_vcf = truth_vcf

        # --- initialize names ---
        # name of the chromosome
        self.chromosome_name = chromosome_name

    def generate_summary(self, start_position, end_position, downsample_rate, bed_list, thread_id):
        """
        Generate labeled images of a given region of the genome
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :param downsample_rate: End position of the region
        :param realignment_flag: If true then realignment will be performed
        :param bed_list: If true then realignment will be performed
        :return:
        """
        if not self.use_hp_info:
            alignment_summarizer = AlignmentSummarizer(self.bam_handler,
                                                       self.fasta_handler,
                                                       self.chromosome_name,
                                                       start_position,
                                                       end_position)

            candidate_images = alignment_summarizer.create_summary(self.truth_vcf,
                                                                   self.train_mode,
                                                                   downsample_rate,
                                                                   bed_list,
                                                                   thread_id,
                                                                   self.random_draw_probability)

            return candidate_images
        else:
            alignment_summarizer_hp = AlignmentSummarizerHP(self.bam_handler,
                                                            self.fasta_handler,
                                                            self.chromosome_name,
                                                            start_position,
                                                            end_position)

            region_summary = alignment_summarizer_hp.create_summary(self.truth_vcf,
                                                                    self.train_mode,
                                                                    downsample_rate,
                                                                    bed_list,
                                                                    self.random_draw_probability)

            if region_summary is not None:
                images_hp1, images_hp2, labels_hp1, labels_hp2, positions, index, image_chunk_ids = region_summary
            else:
                images_hp1, images_hp2, labels_hp1, labels_hp2, positions, index, image_chunk_ids = [], [], [], [], [], [], []

            return images_hp1, images_hp2, labels_hp1, labels_hp2, positions, index, image_chunk_ids



class ImageGenerationUtils:
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
            fasta_handler = PEPPER_VARIANT.FASTA_handler(ref_file)
            bam_handler = PEPPER_VARIANT.BAM_handler(bam_file)
            bam_contigs = bam_handler.get_chromosome_sequence_names()
            fasta_contigs = fasta_handler.get_chromosome_names()
            common_contigs = list(set(fasta_contigs) & set(bam_contigs))
            common_contigs = list(set(common_contigs) - set(EXCLUDED_HUMAN_CONTIGS))

            if len(common_contigs) == 0:
                sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] "
                                 + "ERROR: NO COMMON CONTIGS FOUND BETWEEN THE BAM FILE AND THE FASTA FILE.")
                sys.stderr.flush()
                exit(1)

            common_contigs = sorted(common_contigs, key=ImageGenerationUtils.natural_key)
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: COMMON CONTIGS FOUND: " + str(common_contigs) + "\n")
            sys.stderr.flush()

            for contig_name in common_contigs:
                chromosome_name_list.append((contig_name, None))
        elif chromosome_names:
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
    def generate_image_and_save_to_file(args, all_intervals, total_threads, process_id):
        """
        Method description
        :param args:
        :param all_intervals:
        :param total_threads:
        :param process_id:
        :return:
        """
        thread_prefix = "[THREAD " + "{:02d}".format(process_id) + "]"

        output_path, bam_file, draft_file, use_hp_info, truth_vcf, train_mode, downsample_rate, bed_list, random_draw_probability = args

        timestr = time.strftime("%m%d%Y_%H%M%S")
        file_name = output_path + "pepper_variants_images_thread_" + str(process_id) + "_" + str(timestr)

        if use_hp_info:
            file_name = file_name + "_" + "hp"

        file_name = file_name + ".pkl"

        intervals = [r for i, r in enumerate(all_intervals) if i % total_threads == process_id]

        # initial notification
        if process_id == 0:
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: "
                             + "STARTING PROCESS: " + str(process_id)
                             + " FOR " + str(len(intervals)) + " INTERVALS\n")
        sys.stderr.flush()

        start_time = time.time()
        all_candidates = []
        for counter, interval in enumerate(intervals):
            chr_name, _start, _end = interval

            image_generator = ImageGenerator(chromosome_name=chr_name,
                                             bam_file_path=bam_file,
                                             draft_file_path=draft_file,
                                             use_hp_info=use_hp_info,
                                             truth_vcf=truth_vcf,
                                             train_mode=train_mode,
                                             random_draw_probability=random_draw_probability)

            if not use_hp_info:
                candidates = image_generator.generate_summary(_start, _end, downsample_rate, bed_list, process_id)

                if candidates is not None:
                    # load the previously written objects
                    # if os.path.exists(file_name):
                    #     with gzip.open(file_name, 'rb') as rfp:
                    #         previous_candidates = pickle.load(rfp)
                    #         candidates.extend(previous_candidates)
                    all_candidates.extend(candidates)

                if counter > 0 and counter % 10 == 0 and process_id == 0:
                    percent_complete = int((100 * counter) / len(intervals))
                    time_now = time.time()
                    mins = int((time_now - start_time) / 60)
                    secs = int((time_now - start_time)) % 60

                    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "]"
                                     + " INFO: " + str(thread_prefix) + " " + str(counter) + "/" + str(len(intervals))
                                     + " COMPLETE (" + str(percent_complete) + "%)"
                                     + " [ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec]\n")
                    sys.stderr.flush()

        pickle_output = open(file_name, 'wb')
        pickle.dump(all_candidates, pickle_output)
        pickle_output.close()

        return process_id

    @staticmethod
    def generate_images(chr_list,
                        bam_file,
                        draft_file,
                        use_hp_info,
                        truth_vcf,
                        output_path,
                        total_processes,
                        train_mode,
                        downsample_rate,
                        bed_list,
                        random_draw_probability):
        """
        Description of this method
        :param chr_list:
        :param bam_file:
        :param draft_file:
        :param use_hp_info:
        :param truth_vcf:
        :param output_path:
        :param total_processes:
        :param train_mode:
        :param downsample_rate:
        :param bed_list:
        :return:
        """

        max_size = 100000

        start_time = time.time()
        fasta_handler = PEPPER_VARIANT.FASTA_handler(draft_file)

        all_intervals = []
        total_bases = 0
        # first calculate all the intervals that we need to process
        for chr_name, region in chr_list:
            # contig update message
            if not region:
                interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(chr_name) - 1)
            else:
                interval_start, interval_end = tuple(region)
                interval_start = max(0, interval_start)
                interval_end = min(interval_end, fasta_handler.get_chromosome_sequence_length(chr_name) - 1)

            interval_size = interval_end - interval_start
            if train_mode and interval_size < ImageSizeOptions.MIN_SEQUENCE_LENGTH:
                continue

            # this is the interval size each of the process is going to get which is 10^6
            # I will split this into 10^4 size inside the worker process
            for pos in range(interval_start, interval_end, max_size):
                pos_start = max(interval_start, pos - ImageSizeOptions.MIN_IMAGE_OVERLAP)
                pos_end = min(interval_end, pos + max_size + ImageSizeOptions.MIN_IMAGE_OVERLAP)

                inv_size = pos_end - pos_start
                if train_mode and inv_size < ImageSizeOptions.MIN_SEQUENCE_LENGTH:
                    continue

                all_intervals.append((chr_name, pos_start, pos_end))
                total_bases += inv_size

        # all intervals calculated now
        # contig update message
        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] "
                         + "INFO: TOTAL CONTIGS: " + str(len(chr_list))
                         + " TOTAL INTERVALS: " + str(len(all_intervals))
                         + " TOTAL BASES: " + str(total_bases) + "\n")
        sys.stderr.flush()

        args = (output_path, bam_file, draft_file, use_hp_info, truth_vcf, train_mode, downsample_rate, bed_list, random_draw_probability)
        with concurrent.futures.ProcessPoolExecutor(max_workers=total_processes) as executor:
            futures = [executor.submit(ImageGenerationUtils.generate_image_and_save_to_file, args, all_intervals, total_processes, process_id)
                       for process_id in range(0, total_processes)]

            for fut in concurrent.futures.as_completed(futures):
                if fut.exception() is None:
                    # get the results
                    process_id = fut.result()
                    if process_id == 0:
                        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: THREAD "
                                         + str(process_id) + " FINISHED SUCCESSFULLY.\n")
                else:
                    sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
                fut._result = None  # python issue 27144

        end_time = time.time()
        mins = int((end_time - start_time) / 60)
        secs = int((end_time - start_time)) % 60
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED IMAGE GENERATION\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL ELAPSED TIME FOR GENERATING IMAGES: " + str(mins) + " Min " + str(secs) + " Sec\n")
