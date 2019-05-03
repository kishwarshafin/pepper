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
    if s.lower() not in {'false', 'true', '1', 't'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true' or s.lower() == 't' or s.lower() == '1'


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file containing mapping between reads and the draft assembly."
    )
    parser.add_argument(
        "--draft",
        type=str,
        required=True,
        help="FASTA file containing the draft assembly."
    )
    parser.add_argument(
        "--truth_bam",
        type=str,
        default=None,
        help="BAM file containing mapping of true assembly to the draft assembly."
    )
    parser.add_argument(
        "--chromosome_name",
        type=str,
        help="Region in [chr_name:start-end] format"
    )
    parser.add_argument(
        "--train_mode",
        type=boolean_string,
        default=False,
        help="If true then labeled images will be generated."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="candidate_finder_output/",
        help="Path to output directory, if it doesn't exist it will be created."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=5,
        help="Number of threads to use. Default is 5."
    )
    parser.add_argument(
        "--thread_id",
        type=int,
        required=False,
        help="Thread ID."
    )
    parser.add_argument(
        "--downsample_rate",
        type=float,
        required=False,
        default=1.0,
        help="Downsample reads by this margin."
    )
    FLAGS, unparsed = parser.parse_known_args()
    chr_list = UserInterfaceSupport.get_chromosome_list(FLAGS.chromosome_name, FLAGS.draft)

    if FLAGS.train_mode:
        if not FLAGS.truth_bam:
            raise Exception(TextColor.RED + "ERROR: TRAIN MODE REQUIRES --truth_bam TO BE SET.\n" + TextColor.END)

    output_dir = UserInterfaceSupport.handle_output_directory(os.path.abspath(FLAGS.output_dir))

    UserInterfaceSupport.chromosome_level_parallelization(chr_list,
                                                          FLAGS.bam,
                                                          FLAGS.draft,
                                                          FLAGS.truth_bam,
                                                          output_dir,
                                                          FLAGS.threads,
                                                          FLAGS.thread_id,
                                                          FLAGS.train_mode,
                                                          FLAGS.downsample_rate)
