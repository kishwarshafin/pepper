import argparse
import sys
from datetime import datetime
from pepper.version import __version__
from pepper_variant.modules.argparse.CallVariantsArguments import add_call_variant_arguments
from pepper_variant.modules.argparse.MakeImagesArguments import add_make_images_arguments
from pepper_variant.modules.argparse.RunInferenceArguments import add_run_inference_arguments
from pepper_variant.modules.argparse.FindCandidatesArguments import add_find_candidates_arguments
from pepper_variant.modules.argparse.MergeVariantsArguments import add_merge_variants_arguments
from pepper_variant.modules.argparse.SetParameters import set_parameters
from pepper_variant.modules.python.ImageGenerationUI import ImageGenerationUtils
from pepper_variant.modules.python.RunInference import run_inference
from pepper_variant.modules.python.FindCandidates import process_candidates
from pepper_variant.modules.python.CallVariant import call_variant
from pepper_variant.modules.python.MergeVariants import merge_vcf_records


def main():
    """
    Main interface for PEPPER Variant. The submodules supported as of now are these:
    1) Make images
    2) Call consensus
    3) Stitch
    """
    parser = argparse.ArgumentParser(description='PEPPER VARIANT is the variant calling module of PEPPER.',
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

    parser_merge_variants = subparsers.add_parser('merge_variants', help="Merge SNP variants from PEPPER and DeepVariant.")
    add_merge_variants_arguments(parser_merge_variants)

    options, unparsed = parser.parse_known_args()

    # Following parameters are only used during training, so turning them off by default.
    options.train_mode = False
    options.truth_vcf = None
    options.random_draw_probability = 1.0

    if options.sub_command in ['call_variant', 'make_images', 'find_candidates']:
        options = set_parameters(options)

    if options.sub_command == 'call_variant':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CALL VARIANT MODULE SELECTED\n")
        options.dry = False
        call_variant(options)

    elif options.sub_command == 'make_images':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MAKE IMAGE MODULE SELECTED.\n")
        options.image_output_directory = options.output_dir
        ImageGenerationUtils.generate_images(options)

    elif options.sub_command == 'run_inference':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RUN INFERENCE MODULE SELECTED.\n")
        run_inference(options,
                      options.image_dir,
                      options.output_dir)

    elif options.sub_command == 'find_candidates':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FIND CANDIDATE MODULE SELECTED\n")
        process_candidates(options,
                           options.input_dir,
                           options.output_dir)

    elif options.sub_command == 'merge_variants':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MERGE VARIANTS SUBCOMMAND SELECTED\n")
        merge_vcf_records(options)

    elif options.version is True:
        print("PEPPER VERSION: ", __version__)
    else:
        sys.stderr.write("ERROR: NO SUBCOMMAND SELECTED. PLEASE SELECT ONE OF THE AVAIABLE SUB-COMMANDS.\n")
        parser.print_help()


if __name__ == '__main__':
    main()
