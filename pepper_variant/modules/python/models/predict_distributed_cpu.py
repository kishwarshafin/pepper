import sys
import os
import torch
import time
import torch.onnx
import torch.nn as nn
import onnxruntime
from datetime import datetime
import concurrent.futures
from torch.utils.data import DataLoader
import numpy as np
from numpy import argmax
import h5py

from pepper_variant.modules.python.models.dataloader_predict import SequenceDataset
from pepper_variant.modules.python.models.ModelHander import ModelHandler
from pepper_variant.modules.python.Options import ImageSizeOptions, ImageSizeOptionsHP, TrainOptions
from pepper_variant.modules.python.DataStorePredict import DataStore


def get_all_summary_names(input_file):
    summary_names = []
    with h5py.File(input_file, 'r') as hdf5_file:
        if 'summaries' in hdf5_file:
            summary_names = list(hdf5_file['summaries'].keys())

    return summary_names


def chunks(file_list, n):
    for i in range(0, len(file_list), n):
        yield file_list[i:i + n]


def predict(options, input_filepath, file_chunks, output_filepath, threads, thread_id):
    # create output file
    output_filename = output_filepath + "pepper_prediction_" + str(thread_id) + ".hdf"
    prediction_data_file = DataStore(output_filename, mode='w')

    if thread_id == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: SETTING THREADS TO: " + str(threads) + ".\n")
        sys.stderr.flush()

    # session options
    sess_options = onnxruntime.SessionOptions()
    # sess_options.execution_mode = onnxruntime.ExecutionMode.ORT_SEQUENTIAL
    sess_options.inter_op_num_threads = 1
    sess_options.intra_op_num_threads = 1
    sess_options.execution_mode = onnxruntime.ExecutionMode.ORT_SEQUENTIAL
    sess_options.graph_optimization_level = onnxruntime.GraphOptimizationLevel.ORT_ENABLE_ALL

    if options.quantized:
        ort_session = onnxruntime.InferenceSession(options.quantized_model, sess_options=sess_options)
    else:
        ort_session = onnxruntime.InferenceSession(options.model_path + ".onnx", sess_options=sess_options)

    torch.set_num_threads(1)

    if thread_id == 0:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: STARTING INFERENCE." + "\n")
        sys.stderr.flush()

    batch_completed = 0
    for input_file in file_chunks:
        summary_names = get_all_summary_names(input_file)
        if thread_id == 0:
            sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " +
                             "INFO: TOTAL SUMMARIES: " + str(len(summary_names)) + ".\n")
            sys.stderr.flush()
        chunked_summary_names = list(chunks(summary_names, 1))

        for summary_index, summary_names in enumerate(chunked_summary_names):
            # data loader
            input_data = SequenceDataset([input_filepath], input_file, summary_names)
            data_loader = DataLoader(input_data,
                                     batch_size=options.batch_size,
                                     shuffle=False,
                                     num_workers=options.num_workers,
                                     collate_fn=SequenceDataset.my_collate)

            with torch.no_grad():
                for contigs, positions, depths, candidates, candidate_frequencies, images in data_loader:
                    sys.stderr.flush()
                    # run inference on onnx mode, which takes numpy inputs
                    ort_inputs = {ort_session.get_inputs()[0].name: images.cpu().numpy()}
                    # the return value comes as a list
                    output_type = ort_session.run(None, ort_inputs)
                    output_type = output_type[0]

                    prediction_data_file.write_prediction(batch_completed, contigs, positions, depths, candidates, candidate_frequencies, output_type)

                    batch_completed += 1

            if thread_id == 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " +
                                 "INFO: SUMMARY PROCESSED " + str(summary_index + 1) + "/" + str(len(chunked_summary_names)) + ".\n")
                sys.stderr.flush()

    return thread_id


def predict_pytorch(input_filepath, file_chunks, output_filepath, model_path, batch_size, num_workers, threads):
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: SETTING THREADS TO: " + str(threads) + ".\n")
    torch.set_num_threads(threads)
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: INTER OP-THREAD SET TO: " + str(torch.get_num_threads()) + ".\n")
    sys.stderr.flush()

    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(model_path,
                                                    image_features=ImageSizeOptions.IMAGE_HEIGHT,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS,
                                                    num_type_classes=ImageSizeOptions.TOTAL_TYPE_LABELS)
    transducer_model.eval()
    # create output file
    output_filename = output_filepath + "pepper_prediction" + ".hdf"
    prediction_data_file = DataStore(output_filename, mode='w')

    # data loader
    input_data = SequenceDataset(input_filepath)
    data_loader = DataLoader(input_data,
                             batch_size=batch_size,
                             shuffle=False,
                             num_workers=0,
                             collate_fn=SequenceDataset.my_collate)

    transducer_model.eval()

    batch_completed = 0
    total_batches = len(data_loader)

    with torch.no_grad():
        for contigs, positions, depths, candidates, candidate_frequencies, images in data_loader:
            sys.stderr.flush()
            images = images.type(torch.FloatTensor)

            # run inference
            output_type = transducer_model(images, False)

            output_type = output_type.detach().cpu().numpy()
            # output_type = output_type.detach().cpu().numpy()

            prediction_data_file.write_prediction(batch_completed, contigs, positions, depths, candidates, candidate_frequencies, output_type)
            batch_completed += 1

            if batch_completed % 100 == 0 and batch_completed != 0:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] " + "INFO: BATCHES PROCESSED " + str(batch_completed) + "/" + str(total_batches) + ".\n")
                sys.stderr.flush()


def predict_distributed_cpu(options, filepath, file_chunks, output_filepath, total_callers, threads_per_caller):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param options: Options set for prediction
    :param filepath: Path to image files to predict on
    :param file_chunks: Path to chunked files
    :param output_filepath: Path to output directory
    :param total_callers: Number of callers to start
    :param threads_per_caller: Number of threads per caller.
    :return: Prediction dictionary
    """
    # predict_pytorch(filepath, file_chunks[0],  output_filepath, model_path, batch_size, num_workers, threads_per_caller)
    if options.use_hp_info:
        image_features = ImageSizeOptionsHP.IMAGE_HEIGHT
    else:
        image_features = ImageSizeOptions.IMAGE_HEIGHT

    transducer_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_simple_model_for_training(options.model_path,
                                                    image_features=image_features,
                                                    num_classes=ImageSizeOptions.TOTAL_LABELS,
                                                    num_type_classes=ImageSizeOptions.TOTAL_TYPE_LABELS)
    transducer_model.eval()

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL LOADING TO ONNX\n")
    x = torch.zeros(1, ImageSizeOptions.CANDIDATE_WINDOW_SIZE + 1, image_features)

    if not os.path.isfile(options.model_path + ".onnx"):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: SAVING MODEL TO ONNX\n")
        torch.onnx.export(transducer_model, x,
                          options.model_path + ".onnx",
                          opset_version=11,
                          do_constant_folding=True,
                          input_names=['image'],
                          output_names=['output_type'],
                          dynamic_axes={'image': {0: 'batch_size'},
                                        'output_type': {0: 'batch_size'}})

    if options.quantized:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MODEL QUANTIZATION ENABLED. \n")
        from onnxruntime.quantization import quantize_dynamic, QuantType
        quantize_dynamic(options.model_path + ".onnx", output_filepath + "pepper_model.quantized.onnx", weight_type=QuantType.QUInt8)
        options.quantized_model = output_filepath + "pepper_model.quantized.onnx"
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: QUANTIZED MODEL SAVED.\n")
    start_time = time.time()

    # file_chunks = None
    # thread_id = 0
    # predict(filepath, file_chunks, output_filepath, model_path, batch_size, num_workers, threads_per_caller, thread_id)

    with concurrent.futures.ProcessPoolExecutor(max_workers=total_callers) as executor:
        futures = [executor.submit(predict, options, filepath, file_chunks[thread_id], output_filepath, threads_per_caller, thread_id)
                   for thread_id in range(0, total_callers)]

        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                # get the results
                thread_id = fut.result()
                if thread_id == 0:
                    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: THREAD "
                                     + str(thread_id) + " FINISHED SUCCESSFULLY.\n")
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    end_time = time.time()
    mins = int((end_time - start_time) / 60)
    secs = int((end_time - start_time)) % 60
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: FINISHED PREDICTION\n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ELAPSED TIME: " + str(mins) + " Min " + str(secs) + " Sec\n")

