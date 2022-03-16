## PEPPER-Margin-DeepVariant usage and parameters

We have made our pipeline customizable by exposing majority of the internal parameters to the main module. This allows better control over the pipeline, however, parameter tuning requires to have a good understanding of the pipeline. If you are unsure about any of the parameters, please [open a github issue](https://github.com/kishwarshafin/pepper/issues/new) for a discussion.

The parameters allow to:
* Use custom PEPPER/Margin/DeepVariant model and replace the default models provided.
* Control PEPPER's candidate finding behavior and tune for higher precision or sensitivity.
* Customize the pipeline by skipping optional parts of the pipeline.


### Basic command
The basic command to run our variant calling pipeline is:
```bash
docker run --ipc=host \
-v "INPUT_DIR":/input \
-v "OUTPUT_DIR":/output \
-u `id -u $USER`:`id -g $USER` \
kishwars/pepper_deepvariant:r0.8 \
run_pepper_margin_deepvariant call_variant \
-b </input/READS_2_REFERENCE.bam> \
-f </input/REF.fasta> \
-o </output/OUTPUT_DIR/> \
-t <THREADS> \
--ont_r9_guppy5_sup
```

### Help message
To see the full help message, please run the following:
```bash
docker run --ipc=host \
-v "INPUT_DIR":/input \
-v "OUTPUT_DIR":/output \
-u `id -u $USER`:`id -g $USER` \
kishwars/pepper_deepvariant:r0.8 \
run_pepper_margin_deepvariant call_variant --help
```

### Required parameters
| Parameter           | Short | Type    | Description                                                          |
|---------------------|-------|---------|----------------------------------------------------------------------|
| --bam               | -b    | String  | Path to a bam file containing mapping between reads and a reference. |
| --fasta             | -f    | String  | Path to a reference file in FASTA format.                            |
| --output_dir        | -o    | String  | Path to a output directory.                                          |
| --threads           | -t    | Integer | Number of threads to use.                                            |
| --ont_r9_guppy5_sup |       | Boolean | Set to call variants on R9.4.1 Guppy 5+ sup Oxford Nanopore reads.   |
| --ont_r10_q20       |       | Boolean | Set to call variants on R10.4 Q20 Oxford Nanopore reads.             |
| --hifi              |       | Boolean | Set to call variants on PacBio HiFi reads.                           |

### Optional parameters
| Parameter                      | Short | Type    | Description                                                                                                     |
|--------------------------------|-------|---------|-----------------------------------------------------------------------------------------------------------------|
| --region                       | -r    | String  | Region in [contig_name:start-end] format.                                                                       |
| --output_prefix                | -p    | String  | Prefix for output filename. Do not include extension in prefix.                                                 |
| --keep_intermediate_bam_files  | -k    | Boolean | If set then intermediate haplotagged bam files will be kept in the intermediate files directory. Default: False |
| --sample_name                  | -s    | String  | Name of the sample to put in the output VCF.                                                                    |
| --gvcf                         |       | Boolean | If set then a gVCF output will be generated.                                                                    |
| --only_pepper                  |       | Boolean | If set then pipeline will stop after PEPPER.                                                                    |
| --only_haplotag                |       | Boolean | If set then pipeline will stop after Margin haplotag stage.                                                     |
| --phased_output                |       | Boolean | If set then Margin phase will generate a phased VCF and BAM at the end when --phased_output is set.             |
| --skip_final_phased_bam        |       | Boolean | If set with phased output then the final output will not have a bam file when --phased_output is set.           |
| --dry                          |       | Boolean | If set then only the commands will be printed. [Debugging]                                                      |
| --help                         | -h    | Boolean | Show this text and exit.                                                                                        |

### Parameters for PEPPER

| Parameter                                    | Type            | Description                                                                                                                      |
|----------------------------------------------|-----------------|----------------------------------------------------------------------------------------------------------------------------------|
| --pepper_model                               | String          | Path to a custom PEPPER model.                                                                                                   |
| --pepper_quantized                           | Boolean         | If set then use quantization for inference while on CPU inference mode. Speeds up inference. May add non-deterministic behavior. |
| --no_pepper_quantized                        | Boolean         | If set then do not use quantization for inference while on CPU inference mode.                                                   |
| --pepper_downsample_rate                     | Float [0.0-1.0] | Downsample rate of reads while generating images. Default is 1.0                                                                 |
| --pepper_region_size                         | Integer         | Region size in bp used to chunk the genome. Default is 100000.                                                                   |
| --pepper_include_supplementary               | Boolean         | If set then supplementary reads will be used. Default is False.                                                                  |
| --pepper_min_mapq                            | Integer         | Minimum mapping quality for read to be considered valid. Default is 5                                                            |
| --pepper_min_snp_baseq                       | Integer         | Minimum base quality for base to be considered valid for SNP.                                                                    |
| --pepper_min_indel_baseq                     | Integer         | Minimum base quality for base to be considered valid for INDELs.                                                                 |
| --pepper_snp_frequency                       | Float [0.0-1.0] | Minimum SNP frequency for a site to be considered to have a variant.                                                             |
| --pepper_insert_frequency                    | Float [0.0-1.0] | Minimum insert frequency for a site to be considered to have a variant.                                                          |
| --pepper_delete_frequency                    | Float [0.0-1.0] | Minimum delete frequency for a site to be considered to have a variant.                                                          |
| --pepper_min_coverage_threshold              | Integer         | Minimum coverage of a site for finding candidates.                                                                               |
| --pepper_candidate_support_threshold         | Integer         | Minimum number of reads supporting a variant to be considered as a candidate.                                                    |
| --pepper_snp_candidate_frequency_threshold   | Float [0.0-1.0] | Minimum frequency for a SNP candidate to be considered to be a candidate variant.                                                |
| --pepper_indel_candidate_frequency_threshold | Float [0.0-1.0] | Minimum frequency for an INDEL candidate to be considered to be a candidate variant.                                             |
| --pepper_skip_indels                         | Boolean         | If set then INDEL calling will be skipped.                                                                                       |
| --pepper_allowed_multiallelics               | Integer         | Max allowed multiallelic variants per site.                                                                                      |
| --pepper_snp_p_value                         | Float [0.0-1.0] | Predicted value threshold used for a SNP to be considered a candidate.                                                           |
| --pepper_insert_p_value                      | Float [0.0-1.0] | Predicted value threshold used for a insert to be considered a candidate.                                                        |
| --pepper_delete_p_value                      | Float [0.0-1.0] | Predicted value threshold used for a delete to be considered a candidate.                                                        |
| --pepper_snp_p_value_in_lc                   | Float [0.0-1.0] | Predicted value used for a SNP to be considered a candidate in low complexity regions.                                           |
| --pepper_insert_p_value_in_lc                | Float [0.0-1.0] | Predicted value used for an insert to be considered a candidate in low complexity regions.                                       |
| --pepper_delete_p_value_in_lc                | Float [0.0-1.0] | Predicted value used for a delete to be considered a candidate in low complexity regions.                                        |
| --pepper_snp_q_cutoff                        | Integer         | GQ cutoff for a SNP variant to be re-genotyped with DeepVariant.                                                                 |
| --pepper_indel_q_cutoff                      | Integer         | GQ cutoff for an INDEL variant to be re-genotyped with DeepVariant.                                                              |
| --pepper_snp_q_cutoff_in_lc                  | Integer         | GQ cutoff for a SNP variant in low complexity region to be re-genotyped with DeepVariant.                                        |
| --pepper_indel_q_cutoff_in_lc                | Integer         | GQ cutoff for an INDEL variant in low complexity region to be re-genotyped with DeepVariant.                                     |
| --pepper_report_snp_above_freq               | Float [0.0-1.0] | Report all SNPs above frequency for re-genotyping even if the predicted value is low. Set 0 to disable this. [Experimental]      |
| --pepper_report_indel_above_freq             | Float [0.0-1.0] | Report all INDELs above frequency for re-genotyping even if the predicted value is low. Set 0 to disable this. [Experimental]    |

### Parameters for Margin
| Parameter                  | Type   | Description                                  |
|----------------------------|--------|----------------------------------------------|
| --margin_haplotag_model    | String | Path to a custom margin model.               |
| --margin_phase_model       | String | Path to a custom margin model.               |

### Parameters for DeepVariant
| Parameter                     | Type    | Description                                                                                              |
|-------------------------------|---------|----------------------------------------------------------------------------------------------------------|
| --dv_model                    | String  | Path to a custom DeepVariant model. If set, then this model will be used for both SNP and INDEL calling. |
| --dv_model_snp                | String  | Custom DeepVariant model for SNP calling.                                                                |
| --dv_model_indel              | String  | Custom DeepVariant model for INDEL calling.                                                              |
| --dv_alt_aligned_pileup       | String  | alt_align_pileup used for make_examples of --dv_model parameter. [none, rows, diff_channels]             |
| --dv_alt_aligned_pileup_snp   | String  | alt_align_pileup used for make_examples of --dv_model_snp parameter. [none, rows, diff_channels]         |
| --dv_alt_aligned_pileup_indel | String  | alt_align_pileup used for make_examples of --dv_model_indel parameter. [none, rows, diff_channels]       |
| --dv_pileup_image_width       | String  | DeepVariant image width. [none, rows, diff_channels]                                                     |
| --dv_realign_reads            | String  | If true then local read alingment will be performed. [set: true/false]                                   |
| --dv_partition_size           | String  | DeepVariant partition_size used for make_examples.                                                       |
| --dv_min_mapping_quality      | Integer | DeepVariant minimum mapping quality.                                                                     |
| --dv_min_base_quality         | Integer | DeepVariant minimum base quality.                                                                        |
| --dv_vsc_min_fraction_indels  | String  | DeepVariant minimum vsc fraction for indels.                                                             |
| --dv_vsc_min_fraction_snps    | String  | DeepVariant minimum vsc fraction for snps.                                                               |
| --dv_sort_by_haplotypes       | String  | If true then haplotype sorting will be used. [set: true/false]                                           |
| --dv_parse_sam_aux_fields     | String  | If true then auxiliary field parsing is enabled. [set: true/false]                                       |
| --dv_add_hp_channel           | String  | If true then hp channel will be added. [set: true/false]                                                 |
| --dv_use_hp_information       | String  | If true then hp information will be properly used. [set: true/false]                                     |
| --dv_use_multiallelic_mode    | String  | If true multiallelic model will be used during post-processing. [set: true/false]                        |

