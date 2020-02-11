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


def polish_genome_distributed_gpu(image_dir, model_path, batch_size, num_workers, output_dir, device_ids):
    sys.stderr.write(TextColor.GREEN + "INFO: DISTRIBUTED GPU SETUP\n" + TextColor.END)

    if device_ids is None:
        total_gpu_devices = torch.cuda.device_count()
        sys.stderr.write(TextColor.GREEN + "INFO: TOTAL GPU AVAILABLE: " + str(total_gpu_devices) + "\n" + TextColor.END)
        device_ids = [i for i in range(0, total_gpu_devices)]
        total_callers = total_gpu_devices
    else:
        device_ids = [int(i) for i in device_ids.split(',')]
        for device_id in device_ids:
            major_capable, minor_capable = torch.cuda.get_device_capability(device=device_id)
            if major_capable < 0:
                sys.stderr.write(TextColor.RED + "ERROR: GPU DEVICE: " + str(device_id) + " IS NOT CUDA CAPABLE.\n" + TextColor.END)
                sys.stderr.write(TextColor.GREEN + "Try running: $ python3\n"
                                                   ">>> import torch \n"
                                                   ">>> torch.cuda.get_device_capability(device=" + str(device_id) + ")\n" + TextColor.END)
            else:
                sys.stderr.write(TextColor.GREEN + "INFO: CAPABILITY OF GPU#" + str(device_id) +":\t" + str(major_capable)
                                 + "-" + str(minor_capable) + "\n" + TextColor.END)
        total_callers = len(device_ids)

    sys.stderr.write(TextColor.GREEN + "INFO: AVAILABLE GPU DEVICES: " + str(device_ids) + "\n" + TextColor.END)

    if total_callers == 0:
        sys.stderr.write(TextColor.RED + "ERROR: NO GPU AVAILABLE BUT GPU MODE IS SET\n" + TextColor.END)
        exit()

    # chunk the inputs
    input_files = get_file_paths_from_directory(image_dir)

    file_chunks = [[] for i in range(total_callers)]
    for i in range(0, len(input_files)):
        file_chunks[i % total_callers].append(input_files[i])

    total_callers = min(total_callers, len(file_chunks))
    sys.stderr.write(TextColor.GREEN + "INFO: TOTAL THREADS: " + str(total_callers) + "\n" + TextColor.END)
    predict_distributed_gpu(image_dir, file_chunks, output_dir, model_path, batch_size, device_ids, num_workers)
    sys.stderr.flush()
    sys.stderr.write(TextColor.GREEN + "\nINFO: " + TextColor.END + "PREDICTION GENERATED SUCCESSFULLY.\n")


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser(description="2_pepper_call_consensus.py performs inference on images "
                                                 "using a trained model.",
                                     formatter_class=argparse.RawTextHelpFormatter)
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
        "--gpu",
        default=False,
        action='store_true',
        help="If set then PyTorch will use GPUs for inference. CUDA required."
    )
    parser.add_argument(
        "-d",
        "--distributed",
        default=False,
        action='store_true',
        help="Distributed gpu inference. This mode will enable one caller per GPU."
    )
    parser.add_argument(
        "-d_ids",
        "--device_ids",
        type=str,
        required=False,
        default=None,
        help="List of gpu device ids to use for inference. Only used in distributed setting.\n"
             "Example usage: --device_ids 0,1,2 (this will create three callers in id 'cuda:0, cuda:1 and cuda:2'\n"
             "If none then it will use all available devices."
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

    if FLAGS.gpu and FLAGS.distributed:
        """
        DO DISTRIBUTED GPU INFERENCE. THIS MODE WILL ENABLE ONE MODEL PER GPU
        """
        polish_genome_distributed_gpu(FLAGS.image_dir,
                                      FLAGS.model_path,
                                      FLAGS.batch_size,
                                      FLAGS.num_workers,
                                      FLAGS.output_dir,
                                      FLAGS.device_ids)
    else:
        polish_genome(FLAGS.image_dir,
                      FLAGS.model_path,
                      FLAGS.batch_size,
                      FLAGS.num_workers,
                      FLAGS.output_dir,
                      FLAGS.gpu,
                      FLAGS.onnx_off)
