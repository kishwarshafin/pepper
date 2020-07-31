import sys
import torch
import time
from datetime import datetime
from pepper_hp.modules.python.ImageGenerationUI import UserInterfaceSupport
from pepper_hp.modules.python.models.predict_cpu import predict_distributed_cpu
from pepper_hp.modules.python.models.predict_distributed_gpu import predict_distributed_gpu
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


def polish_genome_distributed_gpu(image_dir, model_path, batch_size, threads, num_workers, output_dir, callers_per_gpu, device_ids):
    start_time = time.time()

    if device_ids is None:
        total_gpu_devices = torch.cuda.device_count()
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL GPU AVAILABLE: " + str(total_gpu_devices) + "\n")
        device_ids = [i for i in range(0, total_gpu_devices)]

        multiplied_device_ids = []
        for device_id in device_ids:
            for i in range(0, callers_per_gpu):
                multiplied_device_ids.append(device_id)
        device_ids = multiplied_device_ids

        total_callers = len(device_ids)
    else:
        device_ids = [int(i) for i in device_ids.split(',')]
        device_ids = list(set(device_ids))
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

        multiplied_device_ids = []
        for device_id in device_ids:
            for i in range(0, callers_per_gpu):
                multiplied_device_ids.append(device_id)
        device_ids = multiplied_device_ids

        total_callers = len(device_ids)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: AVAILABLE GPU DEVICES: " + str(list(set(device_ids))) + "\n")
    if total_callers == 0:
        sys.stderr.write("ERROR: NO GPU AVAILABLE BUT GPU MODE IS SET\n")
        exit()

    # chunk the inputs
    input_files = get_file_paths_from_directory(image_dir)

    file_chunks = [[] for i in range(total_callers)]
    for i in range(0, len(input_files)):
        file_chunks[i % total_callers].append(input_files[i])

    file_chunks = [x for x in file_chunks if x]

    total_callers = min(total_callers, len(file_chunks))
    threads_per_caller = max(1, int(threads/total_callers))
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL CALLERS: " + str(total_callers) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: TOTAL THREADS PER CALLER: " + str(threads_per_caller) + "\n")
    predict_distributed_gpu(image_dir, file_chunks, output_dir, model_path, batch_size, total_callers, threads_per_caller, device_ids, num_workers)
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PREDICTION GENERATED SUCCESSFULLY.\n")

    end_time = time.time()
    mins = int((end_time - start_time) / 60)
    secs = int((end_time - start_time)) % 60
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec\n")


def polish_genome_distributed_cpu(image_dir, model_path, batch_size, threads, num_workers, output_dir):
    start_time = time.time()
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: DISTRIBUTED CPU SETUP.\n")

    # chunk the inputs
    input_files = get_file_paths_from_directory(image_dir)

    callers = threads

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

    end_time = time.time()
    mins = int((end_time - start_time) / 60)
    secs = int((end_time - start_time)) % 60
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] TOTAL ELAPSED TIME FOR INFERENCE: " + str(mins) + " Min " + str(secs) + " Sec\n")


def run_inference(image_dir,
                  model,
                  batch_size,
                  num_workers,
                  output_dir,
                  device_ids,
                  callers_per_gpu,
                  gpu_mode,
                  threads):
    output_dir = UserInterfaceSupport.handle_output_directory(output_dir)

    if gpu_mode:
        polish_genome_distributed_gpu(image_dir,
                                      model,
                                      batch_size,
                                      threads,
                                      num_workers,
                                      output_dir,
                                      callers_per_gpu,
                                      device_ids)
    else:
        polish_genome_distributed_cpu(image_dir,
                                      model,
                                      batch_size,
                                      threads,
                                      num_workers,
                                      output_dir)
