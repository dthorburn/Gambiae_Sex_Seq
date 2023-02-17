# Gambiae_Sex_Seq
Developmental isoform and expression analyses using a mixture of long and short read techniques. These pipelines were either developed using [Nextfow](https://www.nextflow.io/) version 22.04.5 (data downloaded 22/08/2022). All pipelines were developed to be used on the Imperial HPC PBS system with the queue limitations in mind. As a result, some of the scripts are BASH as opposed to nextflow due to the support for long queues (72+ hours) being removed.

## Quality Control
The `ReadsQC` pipeline uses [FastP](https://github.com/OpenGene/fastp), [Kraken2](https://ccb.jhu.edu/software/kraken2/), [Fastq-Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/), and [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to assess the quality of the reads and the level of potential contamination. Reads can be extracted for any NCBI taxonomic ID given it was in the database. The kraken2 database was built in 11/2022 using the standard bacteria, virus, human databases with an additional 11 mosquito genomes added (GenBank accessions: GCF_000005575.2, GCF_002204515.2, GCF_013141755.1, GCF_013758885.1, GCF_015732765.1, GCF_016920715.1, GCF_017562075.2, GCF_943734665.1, GCF_943734695.1, GCF_943734745.1, GCF_943734845.2). The Fastq-Screen database built upon the standard database and added the genus *Anopheles* data. 

All QC reports will be compiled into a single [MultiQC](https://multiqc.info/) report. For Illumina PE or SE reads this will include all 4 programs, whereas for Nanopore reads only FastP and Kraken2 will be utilised. See process LR1 for more Nanopore QC. Following the report, reads can be extracted using the `KreakenExtract` flag.   

Steps to run this pipeline:
1. Install a recent version of [nextflow](https://github.com/nextflow-io/nextflow).
2. Clone this repository using `git clone` into your working directory.
3. Update the project directory path and add required (and optional) arguments in the `ReadsQC.sh` script, including the read type [`IlluminaSE`, `IlluminaPE`, or `Nanopore`].
4. Run the pipeline using `qsub ReadsQC.sh`.

## Long Read Analyses

### LR1. Nanopore Basecalling and Fast5 QC
This pipeline extracts fast5 files from a tarball and uses [ONT-Guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) to call bases and emits a single fastq with all passed reads. Quality control plots are achieved through the use of [nanoQC](https://github.com/wdecoster/nanoQC) and [pycoQ](https://github.com/a-slide/pycoQC). If length filtering is needed, tools like seqtk can be used to filter by length, for example. 

Steps to run this pipeline:
1. Obtain the ONT-Guppy singularity container from me, or build the GPU version of [ONT-Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revaf_14dec2018/linux-guppy) into a singularity container called `ONT_Guppy_GPU.sif`.
2. Install a recent version of [nextflow](https://github.com/nextflow-io/nextflow).
3. Download a copy of the appropriate reference genome (i.e., the *Anopheles gambiae* reference genome from [vectorbase](https://vectorbase.org/vectorbase/app/record/dataset/DS_2251b21396)). 
4. Clone this repository using `git clone` into your working directory.
5. Update the project directory path, ephemeral mount path, export the nextflow program path, and add required (and optional) arguments in the `ONT_Basecalling_QC.sh` script.
6. Run the pipeline using `qsub ONT_Basecalling_QC.sh`.
