import sys
from datetime import datetime


def set_parameters(options):
    """
    Given a set of options set parameters specific to a sequencing platform.
    :param options: Set of options
    :return: Options with parameters set
    """

    if options.ont:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ONT VARIANT CALLING MODE SELECTED.\n")
        # image generation
        options.min_mapq = 5
        options.min_baseq = 1
        options.snp_frequency = 0.10
        options.insert_frequency = 0.15
        options.delete_frequency = 0.15
        options.min_coverage_threshold = 5
        options.candidate_support_threshold = 2
        options.candidate_frequency_threshold = 0.10
        options.skip_indels = False
        # candidate finding
        options.allowed_multiallelics = 4
        options.snp_p_value = 0.1
        options.insert_p_value = 0.3
        options.delete_p_value = 0.2
    elif options.hifi:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: HiFi VARIANT CALLING MODE SELECTED.\n")
        # image generation
        options.min_mapq = 1
        options.min_baseq = 1
        options.snp_frequency = 0.05
        options.insert_frequency = 0.05
        options.delete_frequency = 0.05
        options.min_coverage_threshold = 2
        options.candidate_support_threshold = 2
        options.candidate_frequency_threshold = 0.05
        options.skip_indels = False
        # candidate finding
        options.allowed_multiallelics = 4
        options.snp_p_value = 0.4
        options.insert_p_value = 0.4
        options.delete_p_value = 0.4
    elif options.clr:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CLR VARIANT CALLING MODE SELECTED.\n")
        # image generation
        options.min_mapq = 5
        options.min_baseq = 1
        options.snp_frequency = 0.10
        options.insert_frequency = 0.20
        options.delete_frequency = 0.20
        options.min_coverage_threshold = 5
        options.candidate_support_threshold = 2
        options.candidate_frequency_threshold = 0.10
        options.skip_indels = True
        # candidate finding
        options.allowed_multiallelics = 2
        options.snp_p_value = 0.3
        options.insert_p_value = 1.0
        options.delete_p_value = 1.0

    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: THRESHOLDS ARE SET TO: \n")
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN MAPQ:\t\t\t\t" + str(options.min_mapq) + "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN BASEQ:\t\t\t\t" + str(options.min_baseq)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN SNP FREQUENCY:\t\t\t" + str(options.snp_frequency)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN INSERT FREQUENCY:\t\t" + str(options.insert_frequency)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN DELETE FREQUENCY:\t\t" + str(options.delete_frequency)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN COVERAGE THRESHOLD:\t\t" + str(options.min_coverage_threshold)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN CANDIDATE SUPPORT:\t\t" + str(options.candidate_support_threshold)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN CANDIDATE FREQUENCY:\t\t" + str(options.candidate_frequency_threshold)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: SKIP INDEL CANDIDATES:\t\t" + str(options.skip_indels)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MAX ALLOWED CANDIDATE IN ONE SITE:\t" + str(options.allowed_multiallelics)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN SNP PREDICTIVE VALUE:\t\t" + str(options.snp_p_value)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN INSERT PREDICTIVE VALUE:\t" + str(options.insert_p_value)+ "\n" )
    sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN DELETE PREDICTIVE VALUE:\t" + str(options.delete_p_value)+ "\n" )

    return options