import time
import os
import re
import sys
import gzip
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
    def __init__(self, chromosome_name, bam_file_path, fasta_file_path):
        """
        Initialize a manager object
        :param chromosome_name: Name of the chromosome
        :param bam_file_path: Path to the BAM file
        :param fasta_file_path: Path to the reference FASTA file
        """
        # --- initialize handlers ---
        # create objects to handle different files and query
        self.bam_handler = PEPPER_VARIANT.BAM_handler(bam_file_path)
        self.fasta_handler = PEPPER_VARIANT.FASTA_handler(fasta_file_path)

        # --- initialize names ---
        # name of the chromosome
        self.chromosome_name = chromosome_name

    def generate_summary(self, options, start_position, end_position, bed_list, thread_id):
        """
        Generate labeled images of a given region of the genome
        :param options: Options for generating images
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :param bed_list: List of regions from bed file. [GIAB high-confidence regions]
        :param thread_id: Process ID.
        :return:
        """
        if not options.use_hp_info:
            alignment_summarizer = AlignmentSummarizer(self.bam_handler,
                                                       self.fasta_handler,
                                                       self.chromosome_name,
                                                       start_position,
                                                       end_position)

            candidate_images = alignment_summarizer.create_summary(options,
                                                                   bed_list,
                                                                   thread_id)

            return candidate_images
        else:
            alignment_summarizer = AlignmentSummarizerHP(self.bam_handler,
                                                         self.fasta_handler,
                                                         self.chromosome_name,
                                                         start_position,
                                                         end_position)

            candidate_images = alignment_summarizer.create_summary(options,
                                                                   bed_list,
                                                                   thread_id)

            return candidate_images


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
    def generate_image_and_save_to_file(options, all_intervals, bed_list, process_id):
        """
        Method description
        :param options: Image generation options.
        :param all_intervals: All intervals.
        :param bed_list: List of intervals from bed file.
        :param process_id: Process id.
        :return:
        """
        thread_prefix = "[THREAD " + "{:02d}".format(process_id) + "]"


        timestr = time.strftime("%m%d%Y_%H%M%S")
        file_name = options.image_output_directory + "pepper_variants_images_thread_" + str(process_id) + "_" + str(timestr)

        if options.use_hp_info:
            file_name = file_name + "_" + "hp"

        file_name = file_name + ".hdf5"

        intervals = [r for i, r in enumerate(all_intervals) if i % options.threads == process_id]

        # initial notification
        if process_id == 0:
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: "
                             + "STARTING PROCESS: " + str(process_id)
                             + " FOR " + str(len(intervals)) + " INTERVALS\n")
        sys.stderr.flush()

        start_time = time.time()
        # print("Starting thread", thread_prefix)
        with DataStore(file_name, 'w') as output_hdf_file:
            for counter, interval in enumerate(intervals):
                chr_name, _start, _end = interval
                image_generator = ImageGenerator(chromosome_name=chr_name,
                                                 bam_file_path=options.bam,
                                                 fasta_file_path=options.fasta)

                candidates = image_generator.generate_summary(options, _start, _end, bed_list, process_id)
                if candidates is not None:
                    all_contig = []
                    all_position = []
                    all_depth = []
                    all_candidates = []
                    all_candidate_frequency = []
                    all_image_matrix = []
                    all_base_label = []
                    all_type_label = []
                    for i, candidate in enumerate(candidates):
                        all_contig.append(candidate.contig)
                        all_position.append(candidate.position)
                        all_depth.append(candidate.depth)
                        all_candidates.append(candidate.candidates)
                        all_candidate_frequency.append(candidate.candidate_frequency)
                        all_image_matrix.append(candidate.image_matrix)
                        all_base_label.append(candidate.base_label)
                        all_type_label.append(candidate.type_label)

                    summary_name = chr_name + "_" + str(_start) + "_" + str(_end)

                    output_hdf_file.write_summary(summary_name,
                                                  all_contig,
                                                  all_position,
                                                  all_depth,
                                                  all_candidates,
                                                  all_candidate_frequency,
                                                  all_image_matrix,
                                                  all_base_label,
                                                  all_type_label,
                                                  options.train_mode)

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

        return process_id

    @staticmethod
    def generate_images(options):
        """
        Generates images.
        :param options: Option for generating images.
        :return:
        """
        chr_list, bed_list = ImageGenerationUtils.get_chromosome_list(options.region, options.fasta, options.bam, region_bed=options.region_bed)
        options.image_output_directory = ImageGenerationUtils.handle_output_directory(os.path.abspath(options.image_output_directory))

        start_time = time.time()
        fasta_handler = PEPPER_VARIANT.FASTA_handler(options.fasta)

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
            if options.train_mode and interval_size < ImageSizeOptions.MIN_SEQUENCE_LENGTH:
                continue

            # this is the interval size each of the process is going to get which is 10^6
            # I will split this into 10^4 size inside the worker process
            for pos in range(interval_start, interval_end, options.region_size):
                pos_start = max(interval_start, pos)
                pos_end = min(interval_end, pos + options.region_size)

                inv_size = pos_end - pos_start
                if options.train_mode and inv_size < ImageSizeOptions.MIN_SEQUENCE_LENGTH:
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

        with concurrent.futures.ProcessPoolExecutor(max_workers=options.threads) as executor:
            futures = [executor.submit(ImageGenerationUtils.generate_image_and_save_to_file, options, all_intervals, bed_list, process_id)
                       for process_id in range(0, options.threads)]

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
