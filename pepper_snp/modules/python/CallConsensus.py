import sys
import torch
from datetime import datetime
from pepper_snp.modules.python.ImageGenerationUI import UserInterfaceSupport
from pepper_snp.modules.python.models.predict import predict
from pepper_snp.modules.python.models.predict_distributed_cpu import predict_distributed_cpu
from pepper_snp.modules.python.models.predict_distributed_gpu import predict_distributed_gpu
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


def polish_genome(image_dir, model_path, batch_size, threads, num_workers, output_dir, gpu_mode):
    sys.stderr.write("INFO: OUTPUT DIRECTORY: " + output_dir + "\n")
    output_filename = output_dir + "pepper_predictions.hdf"
    predict(image_dir, output_filename, model_path, batch_size, threads, num_workers, gpu_mode)
    sys.stderr.write("INFO: PREDICTION GENERATED SUCCESSFULLY.\n")


def polish_genome_distributed_gpu(image_dir, model_path, batch_size, threads, num_workers, output_dir, device_ids):
    if device_ids is None:
        total_gpu_devices = torch.cuda.device_count()
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL GPU AVAILABLE: " + str(total_gpu_devices) + "\n")
        device_ids = [i for i in range(0, total_gpu_devices)]
        total_callers = total_gpu_devices
    else:
        device_ids = [int(i) for i in device_ids.split(',')]
        for device_id in device_ids:
            major_capable, minor_capable = torch.cuda.get_device_capability(device=device_id)
            if major_capable < 0:
                sys.stderr.write("ERROR: GPU DEVICE: " + str(device_id) + " IS NOT CUDA CAPABLE.\n")
                sys.stderr.write("Try running: $ python3\n"
                                 ">>> import torch \n"
                                 ">>> torch.cuda.get_device_capability(device=" + str(device_id) + ")\n")
            else:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CAPABILITY OF GPU#"
                                 + str(device_id) + ":\t" + str(major_capable)
                                 + "-" + str(minor_capable) + "\n")
        total_callers = len(device_ids)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: AVAILABLE GPU DEVICES: " + str(device_ids) + "\n")

    if total_callers == 0:
        sys.stderr.write("ERROR: NO GPU AVAILABLE BUT GPU MODE IS SET\n")
        exit()

    # chunk the inputs
    input_files = get_file_paths_from_directory(image_dir)

    file_chunks = [[] for i in range(threads)]
    for i in range(0, len(input_files)):
        file_chunks[i % threads].append(input_files[i])

    threads = min(threads, len(file_chunks))
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL THREADS: " + str(threads) + "\n")
    print(threads, device_ids)
    exit()
    predict_distributed_gpu(image_dir, file_chunks, output_dir, model_path, batch_size, threads, num_workers)
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PREDICTION GENERATED SUCCESSFULLY.\n")


def polish_genome_distributed_cpu(image_dir, model_path, batch_size, threads, num_workers, output_dir):
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: DISTRIBUTED CPU SETUP.\n")

    # chunk the inputs
    input_files = get_file_paths_from_directory(image_dir)

    callers = min(8, threads)

    file_chunks = [[] for i in range(callers)]
    for i in range(0, len(input_files)):
        file_chunks[i % callers].append(input_files[i])

    file_chunks = [x for x in file_chunks if x]

    callers = min(callers, len(file_chunks))
    threads_per_caller = max(1, int(threads/callers))

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL CALLERS: " + str(callers) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: THREADS PER CALLER: " + str(threads_per_caller) + "\n")
    predict_distributed_cpu(image_dir, file_chunks, output_dir, model_path, batch_size, callers, threads_per_caller, num_workers)
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PREDICTION FINISHED SUCCESSFULLY. \n")


def call_consensus(image_dir,
                   model,
                   batch_size,
                   num_workers,
                   output_dir,
                   device_ids,
                   gpu_mode,
                   distributed,
                   threads):
    output_dir = UserInterfaceSupport.handle_output_directory(output_dir)

    if distributed is False:
        polish_genome(image_dir,
                      model,
                      batch_size,
                      threads,
                      num_workers,
                      output_dir,
                      gpu_mode)
    else:
        if gpu_mode:
            polish_genome_distributed_gpu(image_dir,
                                          model,
                                          batch_size,
                                          threads,
                                          num_workers,
                                          output_dir,
                                          device_ids)
        else:
            polish_genome_distributed_cpu(image_dir,
                                          model,
                                          batch_size,
                                          threads,
                                          num_workers,
                                          output_dir)
