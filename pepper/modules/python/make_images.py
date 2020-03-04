from pepper.modules.python.TextColor import TextColor
from pepper.build import PEPPER
import os
import sys
from pepper.modules.python.ImageGenerationUI import UserInterfaceSupport


def make_images(bam_filepath, fasta_filepath, region, output_dir, threads):
    """
    GENERATE IMAGES WITHOUT ANY LABELS. THIS IS USED BY pepper.py
    :param bam_filepath: Path to the input bam file
    :param fasta_filepath: Path to the input fasta file
    :param region: Specific region of interest
    :param output_dir: Path to the output directory
    :param threads: Number of threads to use
    :return:
    """
    # check the bam file
    if not os.path.isfile(bam_filepath) or not PEPPER.BAM_handler(bam_filepath):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE BAM FILE.\n" + TextColor.END)
        exit(1)

    # check the fasta file
    if not os.path.isfile(fasta_filepath):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE FASTA FILE.\n" + TextColor.END)
        exit(1)
    # check the output directory
    output_dir = UserInterfaceSupport.handle_output_directory(os.path.abspath(output_dir))

    # check number of threads
    if threads <= 0:
        sys.stderr.write(TextColor.RED + "ERROR: THREAD NEEDS TO BE >=0.\n" + TextColor.END)
        exit(1)

    # get the list of contigs
    contig_list = UserInterfaceSupport.get_chromosome_list(region, fasta_filepath)

    # call the parallelization method to generate images in parallel
    UserInterfaceSupport.chromosome_level_parallelization(contig_list,
                                                          bam_filepath,
                                                          fasta_filepath,
                                                          truth_bam=None,
                                                          output_path=output_dir,
                                                          total_threads=threads,
                                                          train_mode=False)


def make_train_images(bam_filepath, truth_bam_filepath, fasta_filepath, region, output_dir, threads):
    """
    GENERATE IMAGES WITHOUT ANY LABELS. THIS IS USED BY pepper.py
    :param bam_filepath: Path to the input bam file.
    :param truth_bam_filepath: Path to the bam where truth is aligned to the assembly.
    :param fasta_filepath: Path to the input fasta file.
    :param region: Specific region of interest.
    :param output_dir: Path to the output directory.
    :param threads: Number of threads to use.
    :return:
    """
    # check the bam file
    if not os.path.isfile(bam_filepath) or not PEPPER.BAM_handler(bam_filepath):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE BAM FILE.\n" + TextColor.END)
        exit(1)

    # check the truth bam file
    if not os.path.isfile(truth_bam_filepath) or not PEPPER.BAM_handler(truth_bam_filepath):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE TRUTH BAM FILE.\n" + TextColor.END)
        exit(1)

    # check the fasta file
    if not os.path.isfile(fasta_filepath):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE FASTA FILE.\n" + TextColor.END)
        exit(1)
    # check the output directory
    output_dir = UserInterfaceSupport.handle_output_directory(os.path.abspath(output_dir))

    # check number of threads
    if threads <= 0:
        sys.stderr.write(TextColor.RED + "ERROR: THREAD NEEDS TO BE >=0.\n" + TextColor.END)
        exit(1)

    # get the list of contigs
    contig_list = UserInterfaceSupport.get_chromosome_list(region, fasta_filepath)

    # call the parallelization method to generate images in parallel
    UserInterfaceSupport.chromosome_level_parallelization(contig_list,
                                                          bam_filepath,
                                                          fasta_filepath,
                                                          truth_bam=truth_bam_filepath,
                                                          output_path=output_dir,
                                                          total_threads=threads,
                                                          train_mode=True)
