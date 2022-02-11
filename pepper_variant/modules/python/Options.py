class ConsensCandidateFinder(object):
    REGION_SAFE_BASES = 100


class ImageSizeOptions(object):
    IMAGE_HEIGHT = 26
    IMAGE_CHANNELS = 1
    CANDIDATE_WINDOW_SIZE = 32

    TOTAL_LABELS = 28
    TOTAL_TYPE_LABELS = 3
    MIN_SEQUENCE_LENGTH = 1000
    decoded_labels = ["HOM-REF", "HET-ALT", "HOM-ALT"]
    decoded_base_labels = ["RR", "RA", "RC", "RT", "RG", "R*", "R#", "AA", "AC", "AT", "AG", "A*", "A#", "CC", "CT", "CG", "C*", "C#", "TT", "TG", "T*", "T#", "GG", "G*", "G#", "**", "*#", "##"]


class ImageSizeOptionsHP(object):
    IMAGE_HEIGHT = 48
    TOTAL_LABELS = 28
    TOTAL_TYPE_LABELS = 3

    CANDIDATE_WINDOW_SIZE = 20
    IMAGE_CHANNELS = 1
    SEQ_LENGTH = 1000
    SEQ_OVERLAP = 50
    LABEL_LENGTH = SEQ_LENGTH

    MIN_SEQUENCE_LENGTH = 1000
    MIN_IMAGE_OVERLAP = 100


class ReadFilterOptions(object):
    MIN_MAPQ = 5
    MIN_BASEQ = 1
    INCLUDE_SUPPLEMENTARY = False


class TruthFilterOptions(object):
    MIN_MAPQ = 60
    MIN_BASEQ = 0
    INCLUDE_SUPPLEMENTARY = True


class CandidateFinderOptions(object):
    ALLELE_FREQ_THRESHOLD = 0.00

    SNP_ALT_FREQ_COEF = 0
    SNP_NON_REF_PROB_COEF = 0.000214
    SNP_ALLELE_WEIGHT_COEF = 0.994386
    SNP_BIAS_TERM = -1.8e-05
    SNP_THRESHOLD = 0.001
    SNP_FREQ_THRESHOLD = 0.10
    SNP_UPPER_FREQ = 0.5

    INSERT_ALT_FREQ_COEF = 0
    INSERT_NON_REF_PROB_COEF = 0.615082
    INSERT_ALLELE_WEIGHT_COEF = 0.523549
    INSERT_BIAS_TERM = 0.015388
    INSERT_THRESHOLD = 0.02
    IN_FREQ_THRESHOLD = 0.20
    IN_FREQ_LOWER_THRESHOLD = 0.10
    IN_UPPER_FREQ = 0.5

    DELETE_ALT_FREQ_COEF = 0
    DELETE_NON_REF_PROB_COEF = 0.450943
    DELETE_ALLELE_WEIGHT_COEF = 0.136223
    DELETE_BIAS_TERM = 1.3e-05
    DELETE_THRESHOLD = 0.01
    DEL_FREQ_LOWER_THRESHOLD = 0.10
    DEL_FREQ_THRESHOLD = 0.20
    DEL_UPPER_FREQ = 0.5

    SAFE_BASES = 20
    ALT_PROB_THRESHOLD = 0.01


class PEPPERVariantCandidateFinderOptions(object):
    MOST_ALLOWED_CANDIDATES_PER_SITE = 2
    SAFE_BASES = 20
    ALT_PROB_THRESHOLD = 0.1


class TrainOptions(object):
    # these two parameters are important, make sure you are sliding in a way that you cover the full sequence length
    # the training loop breaks when current_index + TRAIN_WINDOW > LAST_INDEX. You may lose information if you don't
    # slide correctly
    TRAIN_WINDOW = 100
    WINDOW_JUMP = 50
    GRU_LAYERS = 1
    HIDDEN_SIZE = 256


class AlingerOptions(object):
    # base and map quality
    ALIGNMENT_SAFE_BASES = 20
    MIN_MAP_QUALITY = 20

    MAX_READS_IN_REGION = 5000
    RANDOM_SEED = 2719747673

