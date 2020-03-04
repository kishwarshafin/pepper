from pepper.modules.python.TextColor import TextColor
from pepper.build import PEPPER
import os
import sys
import torch
import time
from pathlib import Path
from pepper.modules.python.make_images import make_images
from pepper.modules.python.call_consensus import call_consensus
from pepper.modules.python.perform_stitch import perform_stitch


def polish(bam_filepath, fasta_filepath, output_path, threads, region,
           model_path, batch_size, gpu_mode, distributed, device_ids,
           num_workers, callers, threads_per_caller):
    """
    Run all the sub-modules to polish an input assembly.
    """
    # check the bam file
    if not os.path.isfile(bam_filepath) or not PEPPER.BAM_handler(bam_filepath):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE BAM FILE.\n" + TextColor.END)
        exit(1)

    # check the fasta file
    if not os.path.isfile(fasta_filepath):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE FASTA FILE.\n" + TextColor.END)
        exit(1)

    # check the model file
    if not os.path.isfile(model_path):
        sys.stderr.write(TextColor.RED + "ERROR: CAN NOT LOCATE MODEL FILE.\n" + TextColor.END)
        exit(1)

    # check number of threads
    if threads <= 0:
        sys.stderr.write(TextColor.RED + "ERROR: THREAD NEEDS TO BE >=0.\n" + TextColor.END)
        exit(1)

    # check batch_size
    if batch_size <= 0:
        sys.stderr.write(TextColor.RED + "ERROR: batch_size NEEDS TO BE >0.\n" + TextColor.END)
        exit(1)

    # check num_workers
    if num_workers < 0:
        sys.stderr.write(TextColor.RED + "ERROR: num_workers NEEDS TO BE >=0.\n" + TextColor.END)
        exit(1)

    threads_per_caller = int(threads / max(1, callers))
    # check number of threads
    if threads_per_caller <= 0:
        sys.stderr.write(TextColor.RED + "ERROR: THREAD PER CALLER NEEDS TO BE >=0.\n" + TextColor.END)
        exit(1)

    # check if gpu inference can be done
    if gpu_mode:
        if not torch.cuda.is_available():
            sys.stderr.write(TextColor.RED + "ERROR: TORCH IS NOT BUILT WITH CUDA.\n" + TextColor.END)
            sys.stderr.write(TextColor.RED + "SEE TORCH CAPABILITY:\n$ python3\n"
                                             ">>> import torch \n"
                                             ">>> torch.cuda.is_available()\n If true then cuda is avilable"
                             + TextColor.END)
            exit(1)

    # check if all devices are available
    if device_ids is not None:
        device_ids = [int(i) for i in device_ids.split(',')]
        for device_id in device_ids:
            major_capable, minor_capable = torch.cuda.get_device_capability(device=device_id)
            if major_capable < 0:
                sys.stderr.write(TextColor.RED + "ERROR: GPU DEVICE: " + str(device_id) + " IS NOT CUDA CAPABLE.\n" + TextColor.END)
                sys.stderr.write(TextColor.GREEN + "Try running: $ python3\n"
                                                   ">>> import torch \n"
                                                   ">>> torch.cuda.get_device_capability(device=" + str(device_id) + ")\n" + TextColor.END)
                exit(1)
            else:
                sys.stderr.write(TextColor.GREEN + "INFO: CAPABILITY OF GPU#" + str(device_id) +":\t" + str(major_capable)
                                 + "-" + str(minor_capable) + "\n" + TextColor.END)

    timestr = time.strftime("%m%d%Y_%H%M%S")

    # run directories
    if output_path[-1] != "/":
        output_parent_directory = str(Path(output_path).resolve().parents[0]) + "/"
    else:
        output_parent_directory = output_path

    image_output_directory = output_parent_directory + "images_" + str(timestr) + "/"
    prediction_output_directory = output_parent_directory + "predictions_" + str(timestr) + "/"

    sys.stderr.write(TextColor.GREEN + "INFO: RUN-ID: " + str(timestr) + "\n" + TextColor.END)
    sys.stderr.write(TextColor.GREEN + "INFO: IMAGE OUTPUT: " + str(image_output_directory) + "\n" + TextColor.END)

    sys.stderr.write(TextColor.GREEN + "STEP 1: GENERATING IMAGES\n" + TextColor.END)
    # call the parallelization method to generate images in parallel
    make_images(bam_filepath,
                fasta_filepath,
                region,
                image_output_directory,
                threads)

    sys.stderr.write(TextColor.GREEN + "STEP 2: RUNNING INFERENCE\n" + TextColor.END)
    sys.stderr.write(TextColor.GREEN + "INFO: PREDICTION OUTPUT: " + str(prediction_output_directory) + "\n" + TextColor.END)
    call_consensus(image_output_directory,
                   model_path,
                   batch_size,
                   num_workers,
                   prediction_output_directory,
                   device_ids,
                   gpu_mode,
                   distributed,
                   callers,
                   threads_per_caller)

    sys.stderr.write(TextColor.GREEN + "STEP 3: RUNNING STITCH\n" + TextColor.END)
    sys.stderr.write(TextColor.GREEN + "INFO: PREDICTION OUTPUT: " + str(output_parent_directory) + "\n" + TextColor.END)
    perform_stitch(prediction_output_directory,
                   output_path,
                   threads)
