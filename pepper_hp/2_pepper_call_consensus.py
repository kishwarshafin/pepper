import argparse
import sys
from modules.python.TextColor import TextColor
from modules.python.ImageGenerationUI import UserInterfaceSupport
from modules.python.models.predict import predict


def polish_genome(csv_file, model_path, batch_size, threads, num_workers, output_dir, gpu_mode):
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "OUTPUT DIRECTORY: " + output_dir + "\n")
    output_filename = output_dir + "pepper_predictions.hdf"
    predict(csv_file, output_filename, model_path, batch_size, threads, num_workers, gpu_mode)
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "PREDICTION GENERATED SUCCESSFULLY.\n")


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser(description="2_pepper_call_consensus.py performs inference on images "
                                                 "using a trained model.")
    parser.add_argument(
        "-i",
        "--image_dir",
        type=str,
        required=True,
        help="Path to directory containing all HDF5 images."
    )
    parser.add_argument(
        "-m",
        "--model_path",
        type=str,
        required=True,
        help="Path to a trained model."
    )
    parser.add_argument(
        "-b",
        "--batch_size",
        type=int,
        required=False,
        default=128,
        help="Batch size for testing, default is 100. Suggested values: 256/512/1024."
    )
    parser.add_argument(
        "-w",
        "--num_workers",
        type=int,
        required=False,
        default=4,
        help="Number of workers for loading images. Default is 4."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number threads for pytorch."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=False,
        default='output',
        help="Path to the output directory."
    )
    parser.add_argument(
        "-g",
        "--gpu_mode",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs for inference. CUDA required."
    )
    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir = UserInterfaceSupport.handle_output_directory(FLAGS.output_dir)
    polish_genome(FLAGS.image_dir,
                  FLAGS.model_path,
                  FLAGS.batch_size,
                  FLAGS.threads,
                  FLAGS.num_workers,
                  FLAGS.output_dir,
                  FLAGS.gpu_mode)

