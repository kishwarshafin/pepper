import os
from pepper_variant.modules.python.ImageGenerationUI import ImageGenerationUtils


def make_images(bam, fasta, truth_bam_h1, truth_bam_h2, region, region_bed, output_dir, threads, downsample_rate, train_mode):
    chr_list, bed_list = ImageGenerationUtils.get_chromosome_list(region, fasta, bam, region_bed=region_bed)
    output_dir = ImageGenerationUtils.handle_output_directory(os.path.abspath(output_dir))

    ImageGenerationUtils.generate_images(chr_list,
                                         bam,
                                         fasta,
                                         truth_bam_h1=truth_bam_h1,
                                         truth_bam_h2=truth_bam_h2,
                                         output_path=output_dir,
                                         total_processes=threads,
                                         train_mode=train_mode,
                                         downsample_rate=downsample_rate,
                                         bed_list=bed_list)
