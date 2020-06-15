import os
import sys
import torch
import time
from datetime import datetime
from pepper_snp.modules.python.ImageGenerationUI import UserInterfaceSupport
from pepper_snp.modules.python.MakeImages import make_images
from pepper_snp.modules.python.RunInference import run_inference
from pepper_snp.modules.python.FindSNPCandidates import find_candidates
from pepper_snp.build import PEPPER_SNP


def call_variant(bam_filepath,
                 fasta_filepath,
                 output_dir,
                 threads,
                 region,
                 model_path,
                 batch_size,
                 gpu_mode,
                 callers_per_gpu,
                 distributed,
                 device_ids,
                 num_workers,
                 sample_name,
                 probability_threshold):
    """
    Run all the sub-modules to polish an input assembly.
    """
    start_time = time.time()
    # check the bam file
    if not os.path.isfile(bam_filepath) or not PEPPER_SNP.BAM_handler(bam_filepath):
        sys.stderr.write("ERROR: CAN NOT LOCATE BAM FILE.\n")
        exit(1)

    # check the fasta file
    if not os.path.isfile(fasta_filepath):
        sys.stderr.write("ERROR: CAN NOT LOCATE FASTA FILE.\n")
        exit(1)

    # check the model file
    if not os.path.isfile(model_path):
        sys.stderr.write("ERROR: CAN NOT LOCATE MODEL FILE.\n")
        exit(1)

    # check number of threads
    if threads <= 0:
        sys.stderr.write("ERROR: THREAD NEEDS TO BE >=0.\n")
        exit(1)

    # check batch_size
    if batch_size <= 0:
        sys.stderr.write("ERROR: batch_size NEEDS TO BE >0.\n")
        exit(1)

    # check num_workers
    if num_workers < 0:
        sys.stderr.write("ERROR: num_workers NEEDS TO BE >=0.\n")
        exit(1)

    # check if gpu inference can be done
    if gpu_mode:
        if not torch.cuda.is_available():
            sys.stderr.write("ERROR: TORCH IS NOT BUILT WITH CUDA.\n")
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
                sys.stderr.write("ERROR: GPU DEVICE: " + str(device_id) + " IS NOT CUDA CAPABLE.\n")
                sys.stderr.write("Try running: $ python3\n"
                                 ">>> import torch \n"
                                 ">>> torch.cuda.get_device_capability(device=" + str(device_id) + ")\n")
                exit(1)
            else:
                sys.stderr.write("INFO: CAPABILITY OF GPU#" + str(device_id) +":\t" + str(major_capable)
                                 + "-" + str(minor_capable) + "\n")

    timestr = time.strftime("%m%d%Y_%H%M%S")

    output_dir = UserInterfaceSupport.handle_output_directory(output_dir)

    image_output_directory = output_dir + "images_" + str(timestr) + "/"
    prediction_output_directory = output_dir + "predictions_" + str(timestr) + "/"

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: RUN-ID: " + str(timestr) + "\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: IMAGE OUTPUT: " + str(image_output_directory) + "\n")

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] STEP 1: GENERATING IMAGES\n")
    # call the parallelization method to generate images in parallel
    make_images(bam_filepath,
                fasta_filepath,
                region,
                image_output_directory,
                threads)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] STEP 2: RUNNING INFERENCE\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PREDICTION OUTPUT: " + str(prediction_output_directory) + "\n")
    run_inference(image_output_directory,
                  model_path,
                  batch_size,
                  num_workers,
                  prediction_output_directory,
                  device_ids,
                  callers_per_gpu,
                  gpu_mode,
                  distributed,
                  threads)

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] STEP 3: RUNNING STITCH\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: PREDICTION OUTPUT: " + str(output_dir) + "\n")
    find_candidates(prediction_output_directory,
                    fasta_filepath,
                    output_dir,
                    threads,
                    sample_name,
                    probability_threshold)

    end_time = time.time()
    mins = int((end_time - start_time) / 60)
    secs = int((end_time - start_time)) % 60
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] TOTAL ELAPSED TIME FOR VARIANT CALLING: " + str(mins) + " Min " + str(secs) + " Sec\n")
