from pepper.build import PEPPER
import os
import sys
import torch
import time
from pathlib import Path
from datetime import datetime
from pepper.modules.python.ImageGenerationUI import UserInterfaceSupport
from pepper.modules.python.make_images import make_images
from pepper.modules.python.call_consensus import call_consensus
from pepper.modules.python.perform_stitch import perform_stitch


def polish(bam_filepath, fasta_filepath, output_path, threads, region,
           model_path, batch_size, gpu_mode, device_ids, num_workers):
    """
    Run all the sub-modules to polish an input assembly.
    """
    # check the bam file
    if not os.path.isfile(bam_filepath) or not PEPPER.BAM_handler(bam_filepath):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: CAN NOT LOCATE BAM FILE.\n")
        exit(1)

    # check the fasta file
    if not os.path.isfile(fasta_filepath):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: CAN NOT LOCATE FASTA FILE.\n")
        exit(1)

    # check the model file
    if not os.path.isfile(model_path):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: CAN NOT LOCATE MODEL FILE.\n")
        exit(1)

    # check number of threads
    if threads <= 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: THREAD NEEDS TO BE >=0.\n")
        exit(1)

    # check batch_size
    if batch_size <= 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: batch_size NEEDS TO BE >0.\n")
        exit(1)

    # check num_workers
    if num_workers < 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: num_workers NEEDS TO BE >=0.\n")
        exit(1)

    callers = threads
    threads_per_caller = int(threads / max(1, callers))
    # check number of threads
    if threads_per_caller <= 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: THREAD PER CALLER NEEDS TO BE >=0.\n")
        exit(1)

    # check if gpu inference can be done
    if gpu_mode:
        if not torch.cuda.is_available():
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: TORCH IS NOT BUILT WITH CUDA.\n")
            sys.stderr.write("SEE TORCH CAPABILITY:\n$ python3\n"
                             ">>> import torch \n"
                             ">>> torch.cuda.is_available()\n If true then cuda is avilable")
            exit(1)

    # check if all devices are available
    if device_ids is not None:
        device_ids = [int(i) for i in device_ids.split(',')]
        for device_id in device_ids:
            major_capable, minor_capable = torch.cuda.get_device_capability(device=device_id)
            if major_capable < 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: GPU DEVICE: " + str(device_id) + " IS NOT CUDA CAPABLE.\n")
                sys.stderr.write("Try running: $ python3\n"
                                 ">>> import torch \n"
                                 ">>> torch.cuda.get_device_capability(device=" + str(device_id) + ")\n")
                exit(1)
            else:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CAPABILITY OF GPU#" + str(device_id) +":\t" + str(major_capable)
                                 + "-" + str(minor_capable) + "\n")

    timestr = time.strftime("%m%d%Y_%H%M%S")

    # run directories
    output_dir = UserInterfaceSupport.handle_output_directory(output_path)

    image_output_directory = output_dir + "images_" + str(timestr) + "/"
    prediction_output_directory = output_dir + "predictions_" + str(timestr) + "/"

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RUN-ID: " + str(timestr) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: IMAGE OUTPUT: " + str(image_output_directory) + "\n")

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] STEP 1: GENERATING IMAGES\n")
    sys.stderr.flush()
    # call the parallelization method to generate images in parallel
    make_images(bam_filepath,
                fasta_filepath,
                region,
                image_output_directory,
                threads)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] STEP 2: RUNNING INFERENCE\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PREDICTION OUTPUT: " + str(prediction_output_directory) + "\n")
    sys.stderr.flush()
    call_consensus(image_output_directory,
                   model_path,
                   batch_size,
                   num_workers,
                   prediction_output_directory,
                   device_ids,
                   gpu_mode,
                   threads)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] STEP 3: RUNNING STITCH\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: STITCH OUTPUT: " + str(output_dir) + "\n")
    sys.stderr.flush()
    perform_stitch(prediction_output_directory,
                   output_dir,
                   threads)
