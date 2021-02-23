## P.E.P.P.E.R.
[![Build Status](https://travis-ci.com/kishwarshafin/pepper.svg?branch=master)](https://travis-ci.com/kishwarshafin/pepper)

PEPPER is a genome inference module based on recurrent neural networks that enables long read germline variant calling and assembly polishing in the PEPPER-Margin-DeepVariant pipeline. PEPPER-Margin-DeepVariant produces state-of-the-art variant calling results for both PacBio-HiFi and Oxford Nanopore data. Our pipeline is able to polish nanopore-based diploid assemblies with nanopore or PacBio-HiFi data.  

PEPPER-Margin-DeepVariant pipeline is developed in a collaboration between [UC Santa Cruz genomics institute](https://ucscgenomics.soe.ucsc.edu/) and the [Genomics team in Google Health](https://health.google/health-research/genomics/).

### Case-studies: details results and benchmarking included
 * [Oxford nanopore variant calling case-study](./docs/ONT_variant_calling.md) on HG002 chr20.
 * [PacBio-HiFi variant calling case-study](./docs/HiFi_variant_calling.md) on HG002 chr20.
 * Nanopore-based [assembly polishing with nanopore data](./docs/ONT_polishing.md).
 * Nanopore-based [assembly polishing with PacBio-HiFi data](./docs/HiFi_polishing.md).


### How to run PEPPER-Margin-DeepVariant
Oxford Nanopore or PacBio-HiFi based **Variant calling**:
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
Nanopore-based **Assembly polishing** with Oxford Nanopore or PacBio-HiFi data:
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

### Why use PEPPER-Margin-DeepVariant
 * **Accuracy:** PEPPER-Margin-DeepVariant won the [precisionFDA truth challenge v2](https://www.biorxiv.org/content/10.1101/2020.11.13.380741v1) for all benchmarking region and difficult to map region in the Oxford Nanopore category. The haplotype-aware pipeline produces best results for both Oxford Nanopore and PacBio-HiFi long reads.
 * **Speed:** PEPPER-Margin-DeepVariant provides a cheaper and faster solution to PacBio-HiFi haplotype-aware variant calling. The pipeline can be accelerated on GPU.
 * **Usage**: PEPPER-Margin-DeepVariant provides solution to diploid germline variant calling and nanopore-based assembly polishing. It produces the most accurate variant calls and nanopore assemblies.
 * **Phased output**: PEPPER-Margin-DeepVariant can produce high-quality phasing of variants without trio information with nanopore and PacBio-HiFi reads.


### Use Docker container




## Acknowledgement
We are thankful to the developers of these packages: </br>
* [htslib & samtools](http://www.htslib.org/)
* [ssw library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
* [hdf5 python (h5py)](https://www.h5py.org/)
* [pytorch](https://pytorch.org/)

## Fun Fact
<img src="https://vignette.wikia.nocookie.net/marveldatabase/images/7/72/Anthony_Stark_%28Earth-616%29_from_Iron_Man_Vol_5_2_002.jpg/revision/latest?cb=20130407031815" alt="guppy235" width="240p"> <br/>

The name "P.E.P.P.E.R." is inspired from an A.I. created by Tony Stark in the  Marvel Comics (Earth-616). PEPPER is named after Tony Stark's then friend and the CEO of Resilient, Pepper Potts.


© 2021 Kishwar Shafin.
