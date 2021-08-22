# Run PEPPER-Margin-DeepVariant


## Oxford Nanopore based **Variant calling**:
```bash
docker run --ipc=host \
-v "INPUT_DIR":/input \
-v "OUTPUT_DIR":/output \
-u `id -u $USER`:`id -g $USER` \
kishwars/pepper_deepvariant:r0.5 \
run_pepper_margin_deepvariant call_variant \
-b /input/READS_2_REFERENCE.bam \
-f /input/REF.fasta \
-o /output/OUTPUT_DIR/ \
-p OUTPUT_PREFIX \
-t <THREADS> \
--gvcf \ # optional: Produces gVCF output
--phased_output \ # optional: Produces phased output
--ont
```
