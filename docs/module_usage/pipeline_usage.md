# Run PEPPER-Margin-DeepVariant


## Oxford Nanopore based **Variant calling**:
```bash
docker run --ipc=host \
-v "INPUT_DIR":/input \
-v "OUTPUT_DIR":/output \
-u `id -u $USER`:`id -g $USER` \
kishwars/pepper_deepvariant:r0.4 \
run_pepper_margin_deepvariant call_variant \
-b /input/READS_2_REFERENCE.bam \
-f /input/REF.fasta \
-o /output/OUTPUT_DIR/ \
-p OUTPUT_PREFIX \
-t <THREADS> \
--gvcf \ # optional: Produces gVCF output
--phased_output \ # optional: Produces phased output
--ont # use --ccs for PacBio-HiFi data
```

## Nanopore-based **Assembly polishing** with Oxford Nanopore data:
```bash
docker run --ipc=host \
-v "INPUT_DIR":/input \
-v "OUTPUT_DIR":/output \
-u `id -u $USER`:`id -g $USER` \
pepper_deepvariant:r0.4 \
run_pepper_margin_deepvariant polish_assembly \
-b /input/READS_2_ASSEMBLY.bam \
-f /input/ASSEMBLY.fasta \
-o /output/OUTPUT_DIR/ \
-t <THREADS> \
-s <SAMPLE_NAME> \
--ont # use --ccs for PacBio-HiFi reads aligned to the assembly
```
