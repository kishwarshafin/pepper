class ImageSizeOptions(object):
    IMAGE_HEIGHT = 10
    IMAGE_CHANNELS = 1
    SEQ_LENGTH = 1000
    SEQ_OVERLAP = 50
    LABEL_LENGTH = SEQ_LENGTH

    TOTAL_LABELS = 15
    MIN_SEQUENCE_LENGTH = 1000
    MIN_IMAGE_OVERLAP = 100
    decoded_labels = ['**', 'AA', 'AC', 'AT', 'AG', 'A*', 'CC', 'CT', 'CG', 'C*', 'TT', 'TG', 'T*', 'GG', 'G*']
    class_weights =  [0.1,  0.1,  1.0,  1.0,  1.0,  1.0,  0.1,  1.0,  1.0,  1.0,  0.1,  1.0,  1.0,  0.1,  1.0]


class TrainOptions(object):
    # these two parameters are important, make sure you are sliding in a way that you cover the full sequence length
    # the training loop breaks when current_index + TRAIN_WINDOW > LAST_INDEX. You may lose information if you don't
    # slide correctly
    TRAIN_WINDOW = 100
    WINDOW_JUMP = 50
    GRU_LAYERS = 1
    HIDDEN_SIZE = 128


class CandidateOptions(object):
    # any value higher or equal candidate prob will be considered as a candidate
    CANDIDATE_PROB_THRESHOLD = 0.3


class AlingerOptions(object):
    # base and map quality
    ALIGNMENT_SAFE_BASES = 20
    MIN_MAP_QUALITY = 20

    MAX_READS_IN_REGION = 1500
    RANDOM_SEED = 2719747673

