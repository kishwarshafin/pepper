from __future__ import print_function
import numpy as np
from math import log, ceil
from time import time, ctime
import sys
import logging
from datetime import datetime
"""
The Hyperband class used for hyper-parameter optimization.
Optimized from:
https://github.com/zygmuntz/hyperband
"""


class Hyperband:
    """
    Hyper-parameter optimization algorithm implemented
    This class is optimized from the popular hyperband implementation:
    https://github.com/zygmuntz/hyperband
    """
    def __init__(self, get_params_function, try_params_function, max_iteration, downsample_rate, log_directory,
                 model_directory):
        """
        Initialize hyperband object.
        :param get_params_function: A function to get parameters from sample space.
        :param try_params_function: Try a parameter by training and testing it
        :param max_iteration: Maximum iterations per configuration
        :param downsample_rate: Defines configuration downsampling rate
        :param log_directory: Directory where log is saved
        :param model_directory: Directory where model is saved
        """
        self.get_params = get_params_function
        self.try_params = try_params_function

        self.max_iter = max_iteration  # maximum iterations per configuration
        self.eta = downsample_rate  # defines configuration downsampling rate (default = 3)

        self.logeta = lambda x: log(x) / log(self.eta)
        self.s_max = int(self.logeta(self.max_iter))
        self.B = (self.s_max + 1) * self.max_iter

        self.results = []  # list of dicts
        self.counter = 0
        self.best_loss = np.inf
        self.best_acc = 0
        self.best_counter = -1
        self.best_config = ''
        self.log_file = log_directory+'Hyperband_'+datetime.now().strftime("%Y%m%d_%H%M%S")+'.log'
        self.model_dir = model_directory + 'Hyperband_' + datetime.now().strftime("%Y%m%d_%H%M%S") + "_"

        logging.basicConfig(filename=self.log_file, level=logging.INFO)

    # can be called multiple times
    def run(self, skip_last=0):
        """
        The hyper-parameter optimization algorithm.
        :param skip_last: Skip the last iteration
        :return:
        """
        model_id = 0
        for s in reversed(range(self.s_max + 1)):

            # initial number of configurations
            n = int(ceil(self.B / self.max_iter / (s + 1) * self.eta ** s))

            # initial number of iterations per config
            r = self.max_iter * self.eta ** (-s)

            # n random configurations
            model_configs = [(self.get_params(), False, self.model_dir + 'con_' + str(model_id + i) + '.pkl', 0) for i in range(n)]
            model_id += n
            for i in range((s + 1) - int(skip_last)):  # changed from s + 1

                # Run each of the n configs for <iterations>
                # and keep best (n_configs / eta) configurations

                n_configs = n * self.eta ** (-i)
                n_iterations = int(ceil(r * self.eta ** (i)))

                sys.stderr.write(TextColor.BLUE + "\n*** {} configurations x {:.5f} iterations each"
                                 .format(n_configs, n_iterations) + "\n" + TextColor.END)

                logging.info("\n*** {} configurations x {:.1f} iterations each".format(n_configs, n_iterations))

                val_losses = []
                early_stops = []

                for config_index, config in enumerate(model_configs):
                    self.counter += 1
                    sys.stderr.write(TextColor.BLUE + "{} | {} | lowest loss: {} | accuracy: {} | (run {}) | model {}"
                                     .format(self.counter, ctime(), self.best_loss, self.best_acc, self.best_counter,
                                             self.best_config)
                                     + TextColor.END)
                    logging.info("{} | {} | lowest loss so far: {} | (run {})"
                                 .format(self.counter, ctime(), self.best_loss, self.best_counter))

                    start_time = time()

                    logging.info("Iterations:\t" + str(n_iterations))
                    logging.info("Params:\t" + str(config[0]))
                    params, retrain_model, model_path, prev_ite = config
                    transducer_model, model_optimizer, result = self.try_params(n_iterations, config, model_path)

                    assert (type(result) == dict)
                    assert ('loss' in result)

                    seconds = int(round(time() - start_time))
                    sys.stderr.write(TextColor.BLUE + "\n{} seconds.\n".format(seconds) + TextColor.END)

                    loss = result['loss']
                    val_losses.append(loss)

                    early_stop = result.get('early_stop', False)
                    early_stops.append(early_stop)

                    # keeping track of the best result so far (for display only)
                    # could do it be checking results each time, but hey
                    if loss < self.best_loss:
                        self.best_loss = loss
                        self.best_counter = self.counter
                        self.best_acc = result['accuracy']
                        self.best_config = model_path

                    model_configs[config_index] = (params, True, model_path, n_iterations)

                    result['counter'] = self.counter
                    result['seconds'] = seconds
                    result['params'] = params
                    result['model_path'] = model_path
                    result['iterations'] = n_iterations

                    self.results.append(result)

                # select a number of best configurations for the next loop
                # filter out early stops, if any
                indices = np.argsort(val_losses)
                model_configs = [model_configs[i] for i in indices if not early_stops[i]]
                model_configs = model_configs[0:int(n_configs / self.eta)]

        return self.results
