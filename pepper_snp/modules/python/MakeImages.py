import os
from pepper_snp.modules.python.ImageGenerationUI import UserInterfaceSupport


def make_images(bam_file, draft_file, region, output_path, total_threads):
    output_dir = UserInterfaceSupport.handle_output_directory(os.path.abspath(output_path))

    chr_list = UserInterfaceSupport.get_chromosome_list(region, draft_file, bam_file, region_bed=None)

    UserInterfaceSupport.chromosome_level_parallelization(chr_list=chr_list,
                                                          bam_file=bam_file,
                                                          draft_file=draft_file,
                                                          truth_bam_h1=None,
                                                          truth_bam_h2=None,
                                                          output_path=output_dir,
                                                          total_threads=total_threads,
                                                          realignment_flag=False,
                                                          train_mode=False)


def make_train_images(bam_file, draft_file, truth_bam_h1, truth_bam_h2, region, region_bed, output_path, total_threads):
    output_dir = UserInterfaceSupport.handle_output_directory(os.path.abspath(output_path))

    chr_list = UserInterfaceSupport.get_chromosome_list(region, draft_file, bam_file, region_bed=region_bed)

    UserInterfaceSupport.chromosome_level_parallelization(chr_list=chr_list,
                                                          bam_file=bam_file,
                                                          draft_file=draft_file,
                                                          truth_bam_h1=truth_bam_h1,
                                                          truth_bam_h2=truth_bam_h2,
                                                          output_path=output_dir,
                                                          total_threads=total_threads,
                                                          realignment_flag=False,
                                                          train_mode=True)
