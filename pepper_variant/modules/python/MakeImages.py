import os
from pepper_variant.modules.python.ImageGenerationUI import ImageGenerationUtils


def make_images(bam, fasta, use_hp_info, truth_vcf, region, region_bed, output_dir, threads, downsample_rate, train_mode):
    chr_list, bed_list = ImageGenerationUtils.get_chromosome_list(region, fasta, bam, region_bed=region_bed)
    output_dir = ImageGenerationUtils.handle_output_directory(os.path.abspath(output_dir))

    ImageGenerationUtils.generate_images(chr_list,
                                         bam,
                                         fasta,
                                         use_hp_info,
                                         truth_vcf=truth_vcf,
                                         output_path=output_dir,
                                         total_processes=threads,
                                         train_mode=train_mode,
                                         downsample_rate=downsample_rate,
                                         bed_list=bed_list)
