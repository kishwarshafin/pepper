import argparse
import sys
from datetime import datetime
from pepper.version import __version__
from pepper_variant.modules.argparse.CallVariantsArguments import add_call_variant_arguments
from pepper_variant.modules.argparse.MakeImagesArguments import add_make_images_arguments
from pepper_variant.modules.argparse.RunInferenceArguments import add_run_inference_arguments
from pepper_variant.modules.argparse.FindCandidatesArguments import add_find_candidates_arguments
from pepper_variant.modules.python.MakeImages import make_images
from pepper_variant.modules.python.RunInference import run_inference
from pepper_variant.modules.python.FindCandidates import process_candidates
from pepper_variant.modules.python.CallVariant import call_variant


def main():
    """
    Main interface for PEPPER Variant. The submodules supported as of now are these:
    1) Make images
    2) Call consensus
    3) Stitch
    """
    parser = argparse.ArgumentParser(description='PEPPER VARIANT is the variant calling module of PEPPER. This new module does not depend on haplotyping of the reads.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )

    subparsers = parser.add_subparsers(dest='sub_command')
    subparsers.required = False

    parser_call_variant = subparsers.add_parser('call_variant', help="Run the variant calling pipeline. This will run "
                                                                     "make images-> inference -> find_candidates one after another.\n"
                                                                     "The outputs of each step can be run separately using\n"
                                                                     "the appropriate sub-command.")
    add_call_variant_arguments(parser_call_variant)

    parser_make_images = subparsers.add_parser('make_images', help="Generate images that encode summary statistics "
                                                                   "of reads aligned to an assembly.")
    add_make_images_arguments(parser_make_images)

    parser_run_inference = subparsers.add_parser('run_inference', help="Perform inference on generated images using "
                                                                       "a trained model.")
    add_run_inference_arguments(parser_run_inference)

    parser_find_candidates = subparsers.add_parser('find_candidates', help="Find candidate variants.")
    add_find_candidates_arguments(parser_find_candidates)

    FLAGS, unparsed = parser.parse_known_args()

    if FLAGS.sub_command == 'call_variant':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CALL VARIANT MODULE SELECTED\n")
        call_variant(FLAGS.bam,
                     FLAGS.fasta,
                     FLAGS.output_dir,
                     FLAGS.use_hp_info,
                     FLAGS.freq_based,
                     FLAGS.freq,
                     FLAGS.threads,
                     FLAGS.region,
                     FLAGS.model_path,
                     FLAGS.batch_size,
                     FLAGS.gpu,
                     FLAGS.callers_per_gpu,
                     FLAGS.device_ids,
                     FLAGS.num_workers,
                     FLAGS.sample_name,
                     FLAGS.downsample_rate)

    elif FLAGS.sub_command == 'make_images':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MAKE IMAGE MODULE SELECTED.\n")
        make_images(bam=FLAGS.bam,
                    fasta=FLAGS.fasta,
                    use_hp_info=FLAGS.use_hp_info,
                    truth_vcf=None,
                    region=FLAGS.region,
                    region_bed=None,
                    output_dir=FLAGS.output_dir,
                    threads=FLAGS.threads,
                    downsample_rate=FLAGS.downsample_rate,
                    train_mode=False)

    elif FLAGS.sub_command == 'run_inference':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RUN INFERENCE MODULE SELECTED.\n")
        run_inference(FLAGS.image_dir,
                      FLAGS.model_path,
                      FLAGS.use_hp_info,
                      FLAGS.batch_size,
                      FLAGS.num_workers,
                      FLAGS.output_dir,
                      FLAGS.device_ids,
                      FLAGS.callers_per_gpu,
                      FLAGS.gpu,
                      FLAGS.threads,
                      FLAGS.dry)

    elif FLAGS.sub_command == 'find_candidates':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FIND CANDIDATE MODULE SELECTED\n")
        process_candidates(FLAGS.input_dir,
                           FLAGS.fasta,
                           FLAGS.bam,
                           FLAGS.use_hp_info,
                           FLAGS.sample_name,
                           FLAGS.output_dir,
                           FLAGS.threads,
                           FLAGS.freq_based,
                           FLAGS.freq)

    elif FLAGS.version is True:
        print("PEPPER VERSION: ", __version__)
    else:
        sys.stderr.write("ERROR: NO SUBCOMMAND SELECTED. PLEASE SELECT ONE OF THE AVAIABLE SUB-COMMANDS.\n")
        parser.print_help()


if __name__ == '__main__':
    main()
