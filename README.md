# ![nf-core-conva](docs/images/nf-core-conva_logo.png)

**Copy Number Variation analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.03.0--edge-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**nf-core-conva** is a bioinformatics analysis pipeline to infer copy number variation in normal-tumour paired samples using CNVkit tool or cnv_facets tool or both. This pipeline is designed for use with whole genome sequencing data from short-read sequencing platforms like Illumina and Ion Torrent. It takes fastq.gz files (tumour and normal) as an input along with a reference genome in FASTA format and a gene annotation database in RefFlat format [e.g. refFlat.txt for hg38](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/). It also needs a sorted VCF file of common, polymorphic SNPs. For human samples, a good source is the dbSNP file ([common_all_20180418.vcf.gz](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/)) if you are using GRCh38 Human reference genome and its index file in TBI format.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker / Singularity containers making installation trivial and results highly reproducible. It can also be used with Conda packages.

## Pipeline summary

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Quality trimming ([`Cutadapt!`](https://cutadapt.readthedocs.io/en/stable/index.html))
3. Alignment ([`BWA!`](https://github.com/lh3/bwa))
4. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
5. Duplicate read marking ([`picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
6. Infer copy number changes ([`CNVkit`](https://cnvkit.readthedocs.io/en/stable/index.html)) or ([`cnv_facets`](https://github.com/dariober/cnv_facets)) or both
7. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
 
## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [nf-core-docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command in the directory above the nf-core-conva directory:

    ```bash
    nextflow run nf-core-conva -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. It is also highly recommended to use the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) settings to store the images in a central location for future pipeline runs as shown below:
    
    ```bash
    NXF_SINGULARITY_CACHEDIR=/path/to/central/location nextflow run nf-core-conva -profile test,singularity
    ```        

    * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs as shown below:

    ```bash
    NXF_CONDA_CACHEDIR=/path/to/central/location nextflow run nf-core-conva -profile test,conda
    ```

4. Start running your own analysis!


    * Typical command for singularity based CNV analysis (in the directory above the nf-core-conva directory):

        ```bash
         NXF_SINGULARITY_CACHEDIR=/path/to/central/location nextflow run nf-core-conva \
            --input '[/path/to/samplesheet.csv]' \
            --fasta '[/path/to/reference_genome.fa]' \
            --annotationfile '[/path/to/annotaionfile]' \
            --vcf '[/path/to/sorted_vcf_file]' \
            --tbi '[/path/to/vcf_index.tbi]' \
            --tool 'cnvkit or cnvfacets or all' \ # choose any one option
            --outdir '[/path/to/results/folder/]' \
            -profile singularity
        ```
    * Further details about the format of samplesheet.csv can be seen in the [usage](https://github.com/nibscbioinformatics/nf-core-conva/blob/master/docs/usage.md) document.

## Documentation

The nf-core-conva pipeline comes with documentation about the pipeline: [usage](https://github.com/nibscbioinformatics/nf-core-conva/blob/master/docs/usage.md) and [output](https://github.com/nibscbioinformatics/nf-core-conva/blob/master/docs/output.md).

## Credits

nf-core-conva is written by Ravneet Bhuller [@kaurravneet4123](https://github.com/kaurravneet4123).

Many thanks to others who have helped out along the way too, including (but not limited to):
[@MGordon09](https://github.com/MGordon09),
[@MartinFritzsche](https://github.com/MartinFritzsche)

## Citations

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
