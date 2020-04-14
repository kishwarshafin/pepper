import argparse
import os
from modules.python.TextColor import TextColor
from modules.python.ImageGenerationUI import UserInterfaceSupport


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s: string holding boolean value
    :return:
    """
    if s.lower() not in {'false', 'true', '1', 't', '0', 'f'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true' or s.lower() == 't' or s.lower() == '1'


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser(description="1_pepper_make_images.py script generates summary statistics "
                                                 "from the aligment of reads to the draft assembly.")
    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        help="BAM file containing mapping between reads and the draft assembly."
    )
    parser.add_argument(
        "-d",
        "--draft",
        type=str,
        required=True,
        help="FASTA file containing the draft assembly."
    )
    parser.add_argument(
        "-tb",
        "--truth_bam",
        type=str,
        default=None,
        help="BAM file containing mapping of true assembly to the draft assembly."
    )
    parser.add_argument(
        "-hp",
        "--hp_tag",
        type=int,
        required=True,
        default=None,
        help="Haplotype tag to process from the BAM file."
    )
    parser.add_argument(
        "-r",
        "--region",
        type=str,
        default=None,
        help="Region in [chr_name:start-end] format"
    )
    parser.add_argument(
        "-rb",
        "--region_bed",
        default=None,
        type=str,
        help="Region in [chr_name:start-end] format"
    )
    parser.add_argument(
        "-tm",
        "--train_mode",
        type=boolean_string,
        default=False,
        help="If true then labeled images will be generated."
    )
    parser.add_argument(
        "-rf",
        "--realignment_flag",
        type=boolean_string,
        default=False,
        help="If true then local realignment will be performed."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default="candidate_finder_output/",
        help="Path to output directory, if it doesn't exist it will be created."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=5,
        help="Number of threads to use. Default is 5."
    )
    FLAGS, unparsed = parser.parse_known_args()
    chr_list = UserInterfaceSupport.get_chromosome_list(FLAGS.region, FLAGS.draft, FLAGS.region_bed)

    if FLAGS.train_mode:
        if not FLAGS.truth_bam or not FLAGS.truth_bam:
            raise Exception(TextColor.RED + "ERROR: TRAIN MODE REQUIRES --truth_bam TO BE SET.\n" + TextColor.END)

    output_dir = UserInterfaceSupport.handle_output_directory(os.path.abspath(FLAGS.output_dir))

    UserInterfaceSupport.chromosome_level_parallelization(chr_list,
                                                          FLAGS.bam,
                                                          FLAGS.draft,
                                                          FLAGS.truth_bam,
                                                          FLAGS.hp_tag,
                                                          output_dir,
                                                          FLAGS.threads,
                                                          FLAGS.train_mode,
                                                          FLAGS.realignment_flag)
