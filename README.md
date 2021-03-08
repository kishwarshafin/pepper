## P.E.P.P.E.R.
[![Build Status](https://travis-ci.com/kishwarshafin/pepper.svg?branch=master)](https://travis-ci.com/kishwarshafin/pepper)

`PEPPER` is a genome inference module based on recurrent neural networks that enables long-read variant calling and nanopore assembly polishing in the [PEPPER](https://github.com/kishwarshafin/pepper)-[Margin](https://github.com/UCSC-nanopore-cgl/margin)-[DeepVariant](https://github.com/google/deepvariant) pipeline. This pipeline enables nanopore-based variant calling with [DeepVariant](https://github.com/google/deepvariant).

<p align="center">
<img src="./img/PMDV_variant_calling_ONT.png" alt="PEPPER-Margin-DeepVariant Variant Calling Workflow" width="720p"></img>
</p>

---
### Methods and results
A detailed description of `PEPPER-Margin-DeepVariant` methods and results are discussed in the following manuscript:


<details>
<summary><a href="https://www.biorxiv.org/content/10.1101/2021.03.04.433952v1"><b>bioRxiv:</b> Haplotype-aware variant calling enables high accuracy in nanopore long-reads using deep neural networks.</a></summary>
Authors: Kishwar Shafin, Trevor Pesout, Pi-Chuan Chang, Maria Nattestad, Alexey Kolesnikov, Sidharth Goel, <br/> Gunjan Baid, Jordan M Eizenga, Karen H Miga, Paolo Carnevali, Miten Jain, Andrew Carroll, Benedict Paten.
</details>

---
### Quickstart
Please follow the quickstart guides to assess your setup. Please follow case-study documentations for detailed instructions.
* **Docker**: [Oxford Nanopore and PacBio HiFi variant calling quick start](./docs/quickstart/variant_calling_docker_quickstart.md).
* **Singularity**: [Oxford Nanopore and PacBio HiFi variant calling quick start](./docs/quickstart/variant_calling_singularity_quickstart.md).

### Case studies

The variant calling and assembly polishing [pipelines](./docs/module_usage/pipeline_usage.md) can be run on [Docker](https://docs.docker.com/install/linux/docker-ce/ubuntu/) or [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps). The case studies are designed on `chr20` of `HG002` sample.

Please pick the case-study of your pipeline of interest and the associated container runtime Docker or Singularity. The case-studies include input data and benchmarking of the run:

|                       Pipeline                       |                         Docker                         |                               Singularity                               |                     NVIDIA-docker<br>(GPU)                     |
|:----------------------------------------------------:|:------------------------------------------------------:|:-----------------------------------------------------------------------:|:--------------------------------------------------------------:|
|              Nanopore<br>variant calling             |  [Link](./docs/pipeline_docker/ONT_variant_calling.md) |  [Link](./docs/pipeline_singularity/ONT_variant_calling_singularity.md) |  [Link](./docs/pipeline_docker_gpu/ONT_variant_calling_gpu.md) |
|            PacBio HiFi<br>variant calling            | [Link](./docs/pipeline_docker/HiFi_variant_calling.md) | [Link](./docs/pipeline_singularity/HiFi_variant_calling_singularity.md) | [Link](./docs/pipeline_docker_gpu/HiFi_variant_calling_gpu.md) |
|   Nanopore assembly polishing<br>with nanopore data  |     [Link](./docs/pipeline_docker/ONT_polishing.md)    |     [Link](./docs/pipeline_singularity/ONT_polishing_singularity.md)    |     [Link](./docs/pipeline_docker_gpu/ONT_polishing_gpu.md)    |
| Nanopore assembly polishing<br>with PacBio HiFi data |    [Link](./docs/pipeline_docker/HiFi_polishing.md)    |    [Link](./docs/pipeline_singularity/HiFi_polishing_singularity.md)    |    [Link](./docs/pipeline_docker_gpu/HiFi_polishing_gpu.md)    |


#### Use PEPPER or Margin independently
* If you want to run `PEPPER` or `Margin` independent of the pipeline, please follow this [documentation](./docs/module_usage/module_usage.md).
* If you want to install `PEPPER` locally for development, please follow this [documentation](./docs/local_install/install_pepper_locally.md)

#### License
[PEPPER license](./LICENSE), [Margin License](https://github.com/UCSC-nanopore-cgl/margin/blob/master/LICENSE.txt) and [DeepVariant License](https://github.com/google/deepvariant/blob/r1.1/LICENSE) extend to the trained models (PEPPER, Margin and DeepVariant) and container environment (Docker and Singularity).

### Why use PEPPER-Margin-DeepVariant?
 * **Accuracy:** Our pipeline won the [precisionFDA truth challenge v2](https://www.biorxiv.org/content/10.1101/2020.11.13.380741v1) for all benchmarking region and difficult to map region in the Oxford Nanopore category.
 * **Speed:** `PEPPER-Margin-DeepVariant` provides a cheaper and faster solution to PacBio HiFi haplotype-aware variant calling.
 * **Phased output**: `PEPPER-Margin-DeepVariant` can produce high-quality phasing of variants without trio information with nanopore and PacBio HiFi reads.

### Acknowledgement
We are thankful to the developers of these packages:
* [htslib & samtools](http://www.htslib.org/)
* [pytorch](https://pytorch.org/)
* [ONNX](https://onnx.ai/)
* [hdf5 python (h5py)](https://www.h5py.org/)


### Authors
[PEPPER](https://github.com/kishwarshafin/pepper)-[Margin](https://github.com/UCSC-nanopore-cgl/margin)-[DeepVariant](https://github.com/google/deepvariant) pipeline is developed in a collaboration between [UC Santa Cruz genomics institute](https://ucscgenomics.soe.ucsc.edu/) and the [Genomics team in Google Health](https://health.google/health-research/genomics/).


### Fun Fact
<img src="https://vignette.wikia.nocookie.net/marveldatabase/images/7/72/Anthony_Stark_%28Earth-616%29_from_Iron_Man_Vol_5_2_002.jpg/revision/latest?cb=20130407031815" alt="guppy235" width="240p"> <br/>

The name "P.E.P.P.E.R." is inspired from an A.I. created by Tony Stark in the  Marvel Comics (Earth-616).

PEPPER is named after Tony Stark's then friend and the CEO of Resilient, Pepper Potts.
