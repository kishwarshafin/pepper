import argparse
import sys
import torch
from modules.python.TextColor import TextColor
from modules.python.ImageGenerationUI import UserInterfaceSupport
from modules.python.models.predict import predict
from modules.python.models.predict_distributed_gpu import predict_distributed_gpu
from os.path import isfile, join
from os import listdir


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-3:] == 'hdf']
    return file_paths


def polish_genome(csv_file, model_path, batch_size, num_workers, output_dir, gpu_mode, onnx_mode):
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "OUTPUT DIRECTORY: " + output_dir + "\n")
    output_filename = output_dir + "pepper_predictions.hdf"
    predict(csv_file, output_filename, model_path, batch_size, num_workers, gpu_mode, onnx_mode)
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "PREDICTION GENERATED SUCCESSFULLY.\n")


def polish_genome_distributed_gpu(image_dir, model_path, batch_size, num_workers, output_dir):
    sys.stderr.write(TextColor.GREEN + "INFO: DISTRIBUTED GPU SETUP\n" + TextColor.END)
    total_gpu_devices = torch.cuda.device_count()
    sys.stderr.write(TextColor.GREEN + "INFO: TOTAL GPU AVAILABLE: " + str(total_gpu_devices) + "\n" + TextColor.END)
    total_callers = total_gpu_devices

    if threads == 0:
        sys.stderr.write(TextColor.RED + "ERROR: NO GPU AVAILABLE BUT GPU MODE IS SET\n" + TextColor.END)
        exit()

    # chunk the inputs
    input_files = get_file_paths_from_directory(image_dir)

    file_chunks = [[] for i in range(total_callers)]
    for i in range(0, len(input_files)):
        file_chunks[i % total_callers].append(input_files[i])

    threads = min(total_callers, len(file_chunks))
    sys.stderr.write(TextColor.GREEN + "INFO: TOTAL THREADS: " + str(total_callers) + "\n" + TextColor.END)
    predict_distributed_gpu(image_dir, file_chunks, output_dir, model_path, batch_size, total_callers, num_workers)
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
    parser.add_argument(
        "-d",
        "--distributed_on",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs for inference. CUDA required."
    )
    parser.add_argument(
        "-onnx_off",
        "--onnx_off",
        default=True,
        action='store_false',
        help="Turn off cpu acceleration mode (Disabled when GPU is in use)."
    )
    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir = UserInterfaceSupport.handle_output_directory(FLAGS.output_dir)

    if FLAGS.gpu_mode and FLAGS.distributed_on:
        polish_genome_distributed_gpu(FLAGS.image_dir,
                                      FLAGS.model_path,
                                      FLAGS.batch_size,
                                      FLAGS.num_workers,
                                      FLAGS.output_dir)
    else:
        polish_genome(FLAGS.image_dir,
                      FLAGS.model_path,
                      FLAGS.batch_size,
                      FLAGS.num_workers,
                      FLAGS.output_dir,
                      FLAGS.gpu_mode,
                      FLAGS.onnx_off)
