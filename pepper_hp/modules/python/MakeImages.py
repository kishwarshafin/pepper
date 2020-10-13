import os
from pepper_hp.modules.python.ImageGenerationUI import UserInterfaceSupport


def make_images(bam, fasta, region, output_dir, threads, downsample_rate):
    chr_list = UserInterfaceSupport.get_chromosome_list(region, fasta,  bam, region_bed=None)
    output_dir = UserInterfaceSupport.handle_output_directory(os.path.abspath(output_dir))

    UserInterfaceSupport.chromosome_level_parallelization(chr_list,
                                                          bam,
                                                          fasta,
                                                          truth_bam_h1=None,
                                                          truth_bam_h2=None,
                                                          output_path=output_dir,
                                                          total_threads=threads,
                                                          train_mode=False,
                                                          realignment_flag=False,
                                                          downsample_rate=downsample_rate)


def make_train_images(bam, fasta, truth_bam_h1, truth_bam_h2, region, region_bed, output_dir, threads, downsample_rate):
    chr_list = UserInterfaceSupport.get_chromosome_list(region, fasta,  bam, region_bed=region_bed)
    output_dir = UserInterfaceSupport.handle_output_directory(os.path.abspath(output_dir))

    UserInterfaceSupport.chromosome_level_parallelization(chr_list,
                                                          bam,
                                                          fasta,
                                                          truth_bam_h1=truth_bam_h1,
                                                          truth_bam_h2=truth_bam_h2,
                                                          output_path=output_dir,
                                                          total_threads=threads,
                                                          train_mode=True,
                                                          realignment_flag=False,
                                                          downsample_rate=downsample_rate)
