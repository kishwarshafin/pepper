import sys
from datetime import datetime


def set_parameters(options):
    """
    Given a set of options set parameters specific to a sequencing platform.
    :param options: Set of options
    :return: Options with parameters set
    """

    if options.ont_r9_guppy5_sup:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ONT VARIANT CALLING MODE SELECTED.\n")
        # image generation
        if options.sub_command in ['call_variant', 'make_images', 'make_train_images']:
            if options.min_mapq is None:
                options.min_mapq = 5
            if options.min_snp_baseq is None:
                options.min_snp_baseq = 1
            if options.min_indel_baseq is None:
                options.min_indel_baseq = 1
            if options.snp_frequency is None:
                options.snp_frequency = 0.10
            if options.insert_frequency is None:
                options.insert_frequency = 0.15
            if options.delete_frequency is None:
                options.delete_frequency = 0.15
            if options.min_coverage_threshold is None:
                options.min_coverage_threshold = 3
            if options.candidate_support_threshold is None:
                options.candidate_support_threshold = 2
            if options.snp_candidate_frequency_threshold is None:
                options.snp_candidate_frequency_threshold = 0.10
            if options.indel_candidate_frequency_threshold is None:
                options.indel_candidate_frequency_threshold = 0.10
            if not options.skip_indels:
                options.skip_indels = False
        # candidate finding
        if options.sub_command in ['call_variant', 'find_candidates']:
            if options.allowed_multiallelics is None:
                options.allowed_multiallelics = 4
            if options.snp_p_value is None:
                options.snp_p_value = 0.2
            if options.insert_p_value is None:
                options.insert_p_value = 0.25
            if options.delete_p_value is None:
                options.delete_p_value = 0.25
            if options.snp_q_cutoff is None:
                options.snp_q_cutoff = 15
            if options.indel_q_cutoff is None:
                options.indel_q_cutoff = 10
            if options.report_snp_above_freq is None:
                options.report_snp_above_freq = 0.20
            if options.report_indel_above_freq is None:
                options.report_indel_above_freq = 0.30
    if options.ont_r10_q20:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: ONT VARIANT CALLING MODE SELECTED.\n")
        # image generation
        if options.sub_command in ['call_variant', 'make_images', 'make_train_images']:
            if options.min_mapq is None:
                options.min_mapq = 5
            if options.min_snp_baseq is None:
                options.min_snp_baseq = 5
            if options.min_indel_baseq is None:
                options.min_indel_baseq = 5
            if options.snp_frequency is None:
                options.snp_frequency = 0.12
            if options.insert_frequency is None:
                options.insert_frequency = 0.12
            if options.delete_frequency is None:
                options.delete_frequency = 0.12
            if options.min_coverage_threshold is None:
                options.min_coverage_threshold = 3
            if options.candidate_support_threshold is None:
                options.candidate_support_threshold = 2
            if options.snp_candidate_frequency_threshold is None:
                options.snp_candidate_frequency_threshold = 0.10
            if options.indel_candidate_frequency_threshold is None:
                options.indel_candidate_frequency_threshold = 0.12
            if not options.skip_indels:
                options.skip_indels = False
        # candidate finding
        if options.sub_command in ['call_variant', 'find_candidates']:
            if options.allowed_multiallelics is None:
                options.allowed_multiallelics = 3
            if options.snp_p_value is None:
                options.snp_p_value = 0.1
            if options.insert_p_value is None:
                options.insert_p_value = 0.1
            if options.delete_p_value is None:
                options.delete_p_value = 0.1
            if options.snp_q_cutoff is None:
                options.snp_q_cutoff = 5
            if options.indel_q_cutoff is None:
                options.indel_q_cutoff = 10
            if options.report_snp_above_freq is None:
                options.report_snp_above_freq = 0
            if options.report_indel_above_freq is None:
                options.report_indel_above_freq = 0
    elif options.hifi:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: HiFi VARIANT CALLING MODE SELECTED.\n")
        # image generation
        if options.sub_command in ['call_variant', 'make_images', 'make_train_images']:
            if options.min_mapq is None:
                options.min_mapq = 5
            if options.min_snp_baseq is None:
                options.min_snp_baseq = 10
            if options.min_indel_baseq is None:
                options.min_indel_baseq = 10
            if options.snp_frequency is None:
                options.snp_frequency = 0.10
            if options.insert_frequency is None:
                options.insert_frequency = 0.12
            if options.delete_frequency is None:
                options.delete_frequency = 0.12
            if options.min_coverage_threshold is None:
                options.min_coverage_threshold = 2
            if options.candidate_support_threshold is None:
                options.candidate_support_threshold = 2
            if options.snp_candidate_frequency_threshold is None:
                options.snp_candidate_frequency_threshold = 0.10
            if options.indel_candidate_frequency_threshold is None:
                options.indel_candidate_frequency_threshold = 0.12
            if not options.skip_indels:
                options.skip_indels = False
        if options.sub_command in ['call_variant', 'find_candidates']:
            # candidate finding
            if options.allowed_multiallelics is None:
                options.allowed_multiallelics = 4
            if options.snp_p_value is None:
                options.snp_p_value = 0.01
            if options.insert_p_value is None:
                options.insert_p_value = 0.001
            if options.delete_p_value is None:
                options.delete_p_value = 0.01
            if options.snp_q_cutoff is None:
                options.snp_q_cutoff = 20
            if options.indel_q_cutoff is None:
                options.indel_q_cutoff = 50
            if options.report_snp_above_freq is None:
                options.report_snp_above_freq = 0
            if options.report_indel_above_freq is None:
                options.report_indel_above_freq = 0
    # elif options.clr:
    #     sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: CLR VARIANT CALLING MODE SELECTED.\n")
    #     # image generation
    #     if options.sub_command in ['call_variant', 'make_images', 'make_train_images']:
    #         if options.min_mapq is None:
    #             options.min_mapq = 5
    #         if options.min_baseq is None:
    #             options.min_baseq = 1
    #         if options.snp_frequency is None:
    #             options.snp_frequency = 0.10
    #         if options.insert_frequency is None:
    #             options.insert_frequency = 0.12
    #         if options.delete_frequency is None:
    #             options.delete_frequency = 0.12
    #         if options.min_coverage_threshold is None:
    #             options.min_coverage_threshold = 5
    #         if options.candidate_support_threshold is None:
    #             options.candidate_support_threshold = 2
    #         if options.snp_candidate_frequency_threshold is None:
    #             options.snp_candidate_frequency_threshold = 0.10
    #         if options.indel_candidate_frequency_threshold is None:
    #             options.indel_candidate_frequency_threshold = 0.12
    #         if not options.skip_indels:
    #             options.skip_indels = True
    #     # candidate finding
    #     if options.sub_command in ['call_variant', 'find_candidates']:
    #         if options.allowed_multiallelics is None:
    #             options.allowed_multiallelics = 4
    #         if options.snp_p_value is None:
    #             options.snp_p_value = 0.1
    #         if options.insert_p_value is None:
    #             options.insert_p_value = 0.2
    #         if options.delete_p_value is None:
    #             options.delete_p_value = 0.2
    #         if options.snp_q_cutoff is None:
    #             options.snp_q_cutoff = 20
    #         if options.indel_q_cutoff is None:
    #             options.indel_q_cutoff = 20

    if options.sub_command in ['call_variant', 'make_images', 'make_train_images']:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: THRESHOLDS ARE SET TO: \n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN MAPQ:\t\t\t\t" + str(options.min_mapq) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN SNP BASEQ:\t\t\t" + str(options.min_snp_baseq) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN INDEL BASEQ:\t\t\t" + str(options.min_indel_baseq) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN SNP FREQUENCY:\t\t\t" + str(options.snp_frequency) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN INSERT FREQUENCY:\t\t" + str(options.insert_frequency) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN DELETE FREQUENCY:\t\t" + str(options.delete_frequency) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN COVERAGE THRESHOLD:\t\t" + str(options.min_coverage_threshold) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN CANDIDATE SUPPORT:\t\t" + str(options.candidate_support_threshold) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN SNP CANDIDATE FREQUENCY:\t" + str(options.snp_candidate_frequency_threshold) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN INDEL CANDIDATE FREQUENCY:\t" + str(options.indel_candidate_frequency_threshold) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: SKIP INDEL CANDIDATES:\t\t" + str(options.skip_indels) + "\n")
    if options.sub_command in ['call_variant', 'find_candidates']:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MAX ALLOWED CANDIDATE IN ONE SITE:\t" + str(options.allowed_multiallelics) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN SNP PREDICTIVE VALUE:\t\t" + str(options.snp_p_value) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN INSERT PREDICTIVE VALUE:\t" + str(options.insert_p_value) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: MIN DELETE PREDICTIVE VALUE:\t" + str(options.delete_p_value) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: SNP QV CUTOFF FOR RE-GENOTYPING:\t" + str(options.snp_q_cutoff) + "\n")
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: INDEL QV CUTOFF FOR RE-GENOTYPING:\t" + str(options.indel_q_cutoff) + "\n")

    return options