# nf-core/conva: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Preprocessing](#preprocessing)
  * [FastQC](#fastqc) - Raw read QC
  * [Cutadapt](#cutadapt) - Quality trimming
  * [FastQC](#fastqc) - Trimmed read QC
* [Alignment](#alignment)
  * [BWA-mem](#bwa-mem) - Alignment to the reference genome
* [Alignment post-processing](#alignment-post-processing)
  * [SAMtools](#samtools) - Sort and index alignments
  * [picard MarkDuplicates](#picard-markduplicates) - Duplicate read marking
* [Copy number variation detection](#copy-number-vatiation-detection)
  * [CNVkit](#cnvkit) - Infer copy number variation
  * [cnv_facets](#cnvfacets) - Infer copy number variation
* [Quality control](#quality-control)
  * [MultiQC](#multiqc) - Aggregate report describing results and QC from the pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Preprocessing

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `QC/fastqc_raw/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

**NB:** The FastQC plots in this directory are generated relative to the raw, input reads. They may contain regions of low quality.
</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### Cutadapt

<details markdown="1">
<summary>Output files</summary>

* `QC/cutadapt/`
  * `*_fastq.gz`: The trimmed/modified fastq reads. These files are NOT saved in the pipeline, therefore, you will not find them in the directory.
  * `*.cutadapt.log`: Cutadapt log file containing number and percentage of basepairs processed and trimmed.

</details>

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) finds and removes adapter sequences, primers, poly-A tails and other types of low quality sequence from your high-throughput sequencing reads.

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `QC/fastqc_trimmed/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your trimmed fastq files.
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

**NB:** The FastQC plots in this directory are generated relative to the trimmed reads. The regions of low quality have been removed.
</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads.

## Alignment

### BWA-mem

<details markdown="1">
<summary>Output files</summary>

* `Alignment/bwa/index/bwamem2/`
  * `*.{0123,amb,ann,bwt.2bit.64,pac}`: BWA genome index files
* `Alignment/bwa/`
  * `*.bam`: Tumour and normal bam files.

</details>

[BWA-mem](https://github.com/lh3/bwa) is a software package for mapping DNA sequences against a large reference genome, such as the human genome.

## Alignment post-processing

### SAMtools

<details markdown="1">
<summary>Output files</summary>

**NB:** Please note that the SAMtools' sorted and indexed files are NOT saved in the pipeline. Therefore, you won't find them.

* `Alignment/bwa/samtools_stats/`
  * SAMtools `<SAMPLE>.sorted.bam.flagstat`, `<SAMPLE>.sorted.bam.idxstats` and `<SAMPLE>.sorted.bam.stats` files generated from the alignment files.

</details>

### picard MarkDuplicates

<details markdown="1">
<summary>Output files</summary>

* `Alignment/picard/`
  * `<SAMPLE>.markdup.sorted.bam`: Coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM file and so will be saved by default.
  * `<SAMPLE>.markdup.sorted.bam.bai`: BAI index file for coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM index file and so will be saved by default in the given directory.
* `Alignment/picard/samtools_stats/`
  * SAMtools `<SAMPLE>.markdup.sorted.bam.flagstat`, `<SAMPLE>.markdup.sorted.bam.idxstats` and `<SAMPLE>.markdup.sorted.bam.stats` files generated from the duplicate marked alignment files.

* `Alignment/picard/picard_metrics/`
  * `<SAMPLE>.markdup.sorted.MarkDuplicates.metrics.txt`: Metrics file from MarkDuplicates.

</details>

[picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-) locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.

## Copy number variation detection

### CNVkit

<details markdown="1">
<summary>Output files</summary>

* `CNV_results/cnvkit/`
  * `<fastafile>.bed`: File containing the genomic coordinates of all the accessible regions in the given referenc genome.
  * `<fastafile>.target.bed`: File containing the genomic coordinates of the tiled regions used for targeted resequencing or all the regions in case of WGS data. 
  * `<fastafile>.antitarget.bed`: File containing off-target/antitarget regions
  * `<sample>.markdup.sorted.targetcoverage.cnn`: File containing coverage information in the target regions from BAM read depths.
  * `<sample>.markdup.sorted.antitargetcoverage.cnn`: File containing coverage information in the antitarget regions from BAM read depths.
  * `reference.cnn`: File containing reference copy number from the normal samples.
  * `<sample>.markdup.sorted.cnr`: File containing gene copy number ratios.
  * `<sample>.markdup.sorted.cns`: File containing segment's discrete copy number.
  * `<sample>.markdup.sorted.call.cns`: File containing segment's absolute integer copy number.

</details>

[CNVkit](https://cnvkit.readthedocs.io/en/stable/) infers and visualizes copy number from high-throughput DNA sequencing data.
 
### cnv_facets

<details markdown="1">
<summary>Output files</summary>

* `CNV_results/cnvfacets/`
  * `<prefix>.vcf.gz`: VCF file compressed and indexed of copy number variants.
  * `<prefix>.cnv.png`: Summary plot of CNVs across the genome.
  * `<prefix>.cov.pdf`: Histograms of the distribution of read depth (coverage) across all the position in the tumour and normal sample, before and after filtering positions.
  * `<prefix>.spider.pdf`: This is a diagnostic plot to check how well the copy number fits work The estimated segment summaries are plotted as circles where the size of the circle increases with the number of loci in the segment. The expected value for various integer copy number states are drawn as curves for purity ranging from 0 to 0.95. For a good fit, the segment summaries should be close to one of the lines.
  * `<prefix>.csv.gz`: File of nucleotide counts at each SNP in normal and tumour sample.

</details>

[cnv_facets](https://github.com/dariober/cnv_facets) detects somatic copy number variants (CNV) in tumour-normal samples using the [facets](https://github.com/mskcc/facets) package.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `QC/multiqc/`  
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>


[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
