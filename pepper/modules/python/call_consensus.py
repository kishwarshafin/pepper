import sys
import torch
from pepper.modules.python.ImageGenerationUI import UserInterfaceSupport
from pepper.modules.python.TextColor import TextColor
from pepper.modules.python.models.predict import predict
from pepper.modules.python.models.predict_distributed_gpu import predict_distributed_gpu
from pepper.modules.python.models.predict_distributed_cpu import predict_distributed_cpu
from os.path import isfile, join
from os import listdir
import os


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-3:] == 'hdf']
    return file_paths


def polish_genome(csv_file, model_path, batch_size, num_workers, output_dir, gpu_mode):
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "OUTPUT DIRECTORY: " + output_dir + "\n")
    output_filename = output_dir + "pepper_predictions.hdf"
    predict(csv_file, output_filename, model_path, batch_size, num_workers, gpu_mode)
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "PREDICTION GENERATED SUCCESSFULLY.\n")


def polish_genome_distributed_cpu(image_dir, model_path, batch_size, num_workers, output_dir, total_callers, threads):
    sys.stderr.write(TextColor.GREEN + "INFO: DISTRIBUTED CPU SETUP\n" + TextColor.END)

    # chunk the inputs
    input_files = get_file_paths_from_directory(image_dir)

    file_chunks = [[] for i in range(total_callers)]
    for i in range(0, len(input_files)):
        file_chunks[i % total_callers].append(input_files[i])
    file_chunks = [file_chunks[i] for i in range(len(file_chunks)) if len(file_chunks[i]) > 0]

    total_callers = min(total_callers, len(file_chunks))

    sys.stderr.write(TextColor.GREEN + "INFO: SETUP: " + "\n" + TextColor.END)
    sys.stderr.write(TextColor.GREEN + "INFO: TOTAL CALLERS: " + str(total_callers) + "\n" + TextColor.END)
    sys.stderr.write(TextColor.GREEN + "INFO: THREADS PER CALLER: " + str(threads) + "\n" + TextColor.END)
    sys.stderr.write(TextColor.GREEN + "INFO: DATA-LOADER PER CALLER: " + str(num_workers) + "\n" + TextColor.END)
    sys.stderr.flush()
    predict_distributed_cpu(image_dir,
                            file_chunks,
                            output_dir,
                            model_path,
                            batch_size,
                            total_callers,
                            threads,
                            num_workers)
    sys.stderr.flush()
    sys.stderr.write(TextColor.GREEN + "\nINFO: " + TextColor.END + "PREDICTION GENERATED SUCCESSFULLY.\n")


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


def call_consensus(image_dir, model_path, batch_size, num_workers, output_dir, device_ids, gpu, distributed,
                   callers, threads):

    # check the model file
    if not os.path.isfile(model_path):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE MODEL FILE.\n" + TextColor.END)
        exit(1)
    # check the input directory
    if not os.path.isdir(image_dir):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE IMAGE DIRECTORY.\n" + TextColor.END)
        exit(1)

    # check batch_size
    if batch_size <= 0:
        sys.stderr.write(TextColor.RED + "ERROR: batch_size NEEDS TO BE >0.\n" + TextColor.END)
        exit(1)

    # check num_workers
    if num_workers < 0:
        sys.stderr.write(TextColor.RED + "ERROR: num_workers NEEDS TO BE >=0.\n" + TextColor.END)
        exit(1)

    # check number of threads
    if threads <= 0:
        sys.stderr.write(TextColor.RED + "ERROR: THREAD NEEDS TO BE >=0.\n" + TextColor.END)
        exit(1)

    output_dir = UserInterfaceSupport.handle_output_directory(output_dir)

    if gpu:
        """
        DO DISTRIBUTED GPU INFERENCE. THIS MODE WILL ENABLE ONE MODEL PER GPU
        """
        if not torch.cuda.is_available():
            sys.stderr.write(TextColor.RED + "ERROR: TORCH IS NOT BUILT WITH CUDA.\n" + TextColor.END)
            sys.stderr.write(TextColor.RED + "SEE TORCH CAPABILITY:\n$ python3\n"
                                             ">>> import torch \n"
                                             ">>> torch.cuda.is_available()\n If true then cuda is avilable"
                             + TextColor.END)
            exit(1)

        total_gpu_devices = torch.cuda.device_count()

        if distributed and total_gpu_devices > 1:
            # Distributed GPU setup
            polish_genome_distributed_gpu(image_dir,
                                          model_path,
                                          batch_size,
                                          num_workers,
                                          output_dir,
                                          device_ids)
        else:
            # Normal GPU setup
            polish_genome(image_dir,
                          model_path,
                          batch_size,
                          num_workers,
                          output_dir,
                          gpu)
    elif distributed:
        """
        DO DISTRIBUTED CPU INFERENCE. THIS MODE WILL CREATE MULTIPLE CALLERS.
        """
        # distributed CPU setup
        polish_genome_distributed_cpu(image_dir,
                                      model_path,
                                      batch_size,
                                      num_workers,
                                      output_dir,
                                      callers,
                                      threads)
    else:
        # normal CPU setup
        polish_genome(image_dir,
                      model_path,
                      batch_size,
                      num_workers,
                      output_dir,
                      gpu)
