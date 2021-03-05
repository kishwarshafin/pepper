import argparse
import pickle
from collections import defaultdict
import matplotlib.pyplot as plt

def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s: string holding boolean value
    :return:
    """
    if s.lower() not in {'false', 'true', '1', 't', '0', 'f'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true' or s.lower() == 't' or s.lower() == '1'


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pickle_file",
        type=str,
        required=True,
        help="BAM file containing mapping between reads and the draft assembly."
    )
    FLAGS, unparsed = parser.parse_known_args()

    plots_dictionary = defaultdict()
    params_set = set()
    params_epochs = defaultdict(list)

    results = pickle.load(open(FLAGS.pickle_file, "rb"))
    results = sorted(results, key=lambda r: r['loss'])

    for i, result in enumerate(results):
        if result['iterations'] < 100:
            continue
        print(i+1)
        print(result)
        print("Loss:\t\t", result['loss'])
        print("iterations:\t", result['iterations'])
        print("Accuracy:\t", result['accuracy'])
        print("Params:\t\t", result['params'])
        print("Model path:\t", result['model_path'])

    # for i, result in enumerate(results):
    #     # if result['iterations'] < 100:
    #     #         continue
    #     params_set.add((result['params']['lr'], result['params']['l2']))
    #     for epoch, accuracy in result['accuracy_epoch']:
    #         plots_dictionary[(result['params']['lr'], result['params']['l2'], epoch)] = accuracy
    #         params_epochs[(result['params']['lr'], result['params']['l2'])].append(epoch)
    #
    # for param in params_set:
    #     x = params_epochs[param]
    #     if max(x) < 99:
    #         continue
    #
    #     y = list()
    #     for p in x:
    #         y.append(plots_dictionary[(param[0], param[1], p)])
    #     plt.plot(x, y)
    # axes = plt.axes()
    # axes.set_ylim([99.990, 99.995])
    # plt.show()