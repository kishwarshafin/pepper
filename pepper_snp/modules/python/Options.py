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
    # class_weights =  [6095.0, 1.0, 4700.0, 6021.0, 997.0, 16712.0, 1.0, 1018.0, 4109.0, 60284.0, 1.0, 4531.0, 17371.0, 1.0, 62074.0]
    class_weights = [1.0,  0.1,  1.0,  1.0,  1.0,  1.0,  0.1,  1.0,  1.0,  1.0,  0.1,  1.0,  1.0,  0.1,  1.0]


class TrainOptions(object):
    # these two parameters are important, make sure you are sliding in a way that you cover the full sequence length
    # the training loop breaks when current_index + TRAIN_WINDOW > LAST_INDEX. You may lose information if you don't
    # slide correctly
    TRAIN_WINDOW = 100
    WINDOW_JUMP = 50
    GRU_LAYERS = 1
    HIDDEN_SIZE = 128


class Profiles(object):
    # line 274 candidate_finder.cpp
    ONT_PROFILE = 0
    CCS_PROFILE = 1


class CandidateFinderOptions(object):
    MOST_ALLOWED_CANDIDATES_PER_SITE = 2
    SAFE_BASES = 20
    ALT_PROB_THRESHOLD = 0.1


class ReadFilterOptions(object):
    MIN_MAPQ = 1
    MIN_BASEQ = 1
    INCLUDE_SUPPLEMENTARY = False


class AlingerOptions(object):
    # base and map quality
    ALIGNMENT_SAFE_BASES = 20
    MIN_MAP_QUALITY = 1

    MAX_READS_IN_REGION = 5000
    RANDOM_SEED = 2719747673


class ReadOptions(object):
    MIN_MAPPING_QUALITY = 1
    MIN_BASE_QUALITY = 1