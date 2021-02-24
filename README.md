## P.E.P.P.E.R.
[![Build Status](https://travis-ci.com/kishwarshafin/pepper.svg?branch=master)](https://travis-ci.com/kishwarshafin/pepper)

`PEPPER` is a genome inference module based on recurrent neural networks that enables long-read variant calling and nanopore assembly polishing in the `PEPPER-Margin-DeepVariant` pipeline. `PEPPER-Margin-DeepVariant` produces state-of-the-art variant calling results for both PacBio-HiFi and Oxford Nanopore data.

<img src="./img/PMDV_variant_calling_ONT.png" alt="PEPPER-Margin-DeepVariant Variant Calling Workflow" width="720p"></img>

PEPPER-Margin-DeepVariant pipeline is developed in a collaboration between [UC Santa Cruz genomics institute](https://ucscgenomics.soe.ucsc.edu/) and the [Genomics team in Google Health](https://health.google/health-research/genomics/).

### How to run PEPPER-Margin-DeepVariant

`PEPPER-Margin-DeepVariant` [pipelines](./module_usage/pipeline_usage.md) can be run on [Docker](https://docs.docker.com/install/linux/docker-ce/ubuntu/) or [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps).

Please pick the case-study of your pipeline of interest and the associated container runtime Docker or Singularity. The case-studies include input data and benchmarking of the run:

|                       Pipeline                       |                         Docker                         |                               Singularity                               |                     NVIDIA-docker<br>(GPU)                     |
|:----------------------------------------------------:|:------------------------------------------------------:|:-----------------------------------------------------------------------:|:--------------------------------------------------------------:|
|              Nanopore<br>variant calling             |  [Link](./docs/pipeline_docker/ONT_variant_calling.md) |  [Link](./docs/pipeline_singularity/ONT_variant_calling_singularity.md) |  [Link](./docs/pipeline_docker_gpu/ONT_variant_calling_gpu.md) |
|            PacBio-HiFi<br>variant calling            | [Link](./docs/pipeline_docker/HiFi_variant_calling.md) | [Link](./docs/pipeline_singularity/HiFi_variant_calling_singularity.md) | [Link](./docs/pipeline_docker_gpu/HiFi_variant_calling_gpu.md) |
|   Nanopore assembly polishing<br>with nanopore data  |     [Link](./docs/pipeline_docker/ONT_polishing.md)    |     [Link](./docs/pipeline_singularity/ONT_polishing_singularity.md)    |     [Link](./docs/pipeline_docker_gpu/ONT_polishing_gpu.md)    |
| Nanopore assembly polishing<br>with PacBio-HiFi data |    [Link](./docs/pipeline_docker/HiFi_polishing.md)    |    [Link](./docs/pipeline_singularity/HiFi_polishing_singularity.md)    |    [Link](./docs/pipeline_docker_gpu/HiFi_polishing_gpu.md)    |

### Use PEPPER or Margin independently
* If you want to run `PEPPER` or `Margin` independent of the pipeline, please follow this [documentation](./module_usage/module_usage.md).
* If you want to install `PEPPER` locally for development, please follow this [documentation](./local_install/install_pepper_locally.md)


### Why use PEPPER-Margin-DeepVariant
 * **Accuracy:** `PEPPER-Margin-DeepVariant` won the [precisionFDA truth challenge v2](https://www.biorxiv.org/content/10.1101/2020.11.13.380741v1) for all benchmarking region and difficult to map region in the Oxford Nanopore category.
 * **Speed:** `PEPPER-Margin-DeepVariant` provides a cheaper and faster solution to PacBio-HiFi haplotype-aware variant calling.
 * **Phased output**: `PEPPER-Margin-DeepVariant` can produce high-quality phasing of variants without trio information with nanopore and PacBio-HiFi reads.

### Acknowledgement
We are thankful to the developers of these packages: </br>
* [htslib & samtools](http://www.htslib.org/)
* [ssw library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
* [hdf5 python (h5py)](https://www.h5py.org/)
* [pytorch](https://pytorch.org/)

## Fun Fact
<img src="https://vignette.wikia.nocookie.net/marveldatabase/images/7/72/Anthony_Stark_%28Earth-616%29_from_Iron_Man_Vol_5_2_002.jpg/revision/latest?cb=20130407031815" alt="guppy235" width="240p"> <br/>

The name "P.E.P.P.E.R." is inspired from an A.I. created by Tony Stark in the  Marvel Comics (Earth-616).

PEPPER is named after Tony Stark's then friend and the CEO of Resilient, Pepper Potts.
