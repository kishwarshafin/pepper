class ImageSizeOptions(object):
    IMAGE_HEIGHT = 10
    IMAGE_CHANNEL_HEIGHT = 10
    IMAGE_CHANNELS = 1
    SEQ_LENGTH = 1000
    SEQ_OVERLAP = 50
    LABEL_LENGTH = SEQ_LENGTH

    TOTAL_LABELS = 5
    MIN_SEQUENCE_LENGTH = 1000
    MIN_IMAGE_OVERLAP = 100


class ReadFilterOptions(object):
    MIN_MAPQ = 10
    MIN_BASEQ = 0
    INCLUDE_SUPPLEMENTARY = False


class CandidateFinderOptions(object):
    SAFE_BASES = 20
    # ALT_PROB_THRESHOLD = 0.4
    # NON_REF_PROB_THRESHOLD = 0.4
    # ALLELE_FREQ_THRESHOLD = 0.11
    # ALLELE_FREQ_THRESHOLD_LAST_RESORT = 0.1
    # PROB_LAST_RESORT = 0.5

    ALT_PROB_THRESHOLD = 0.1
    NON_REF_PROB_THRESHOLD = 0.1
    ALLELE_FREQ_THRESHOLD = 0.12
    ALLELE_FREQ_THRESHOLD_LAST_RESORT = 0.20
    PROB_LAST_RESORT = 0.1


class TrainOptions(object):
    # these two parameters are important, make sure you are sliding in a way that you cover the full sequence length
    # the training loop breaks when current_index + TRAIN_WINDOW > LAST_INDEX. You may lose information if you don't
    # slide correctly
    TRAIN_WINDOW = 100
    WINDOW_JUMP = 50
    GRU_LAYERS = 1
    HIDDEN_SIZE = 128


class AlingerOptions(object):
    # base and map quality
    ALIGNMENT_SAFE_BASES = 20
    MIN_MAP_QUALITY = 20

    MAX_READS_IN_REGION = 1500
    RANDOM_SEED = 2719747673

