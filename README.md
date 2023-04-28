# Gambiae SexSeq
Developmental isoform and expression analyses using a mixture of long and short read techniques. These pipelines were either developed using [Nextfow](https://www.nextflow.io/) version 20.10.0. All pipelines were developed to be used on the Imperial HPC PBS system with the queue limitations in mind. As a result, some of the scripts are BASH as opposed to nextflow due to the support for long queues (72+ hours) being removed.

## Quality control
The `ReadsQC` pipeline uses [FastP](https://github.com/OpenGene/fastp), [Kraken2](https://ccb.jhu.edu/software/kraken2/), [Fastq-Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/), and [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to assess the quality of the reads and the level of potential contamination. Reads can be extracted for any NCBI taxonomic ID given it was in the database. The kraken2 database was built in 11/2022 using the standard bacteria, virus, human databases with an additional 11 mosquito genomes added (GenBank accessions: GCF_000005575.2, GCF_002204515.2, GCF_013141755.1, GCF_013758885.1, GCF_015732765.1, GCF_016920715.1, GCF_017562075.2, GCF_943734665.1, GCF_943734695.1, GCF_943734745.1, GCF_943734845.2). The Fastq-Screen database built upon the standard database and added the genus *Anopheles* data. 

All QC reports will be compiled into a single [MultiQC](https://multiqc.info/) report. For Illumina PE or SE reads this will include all 4 programs, whereas for Nanopore reads only FastP and Kraken2 will be utilised. See process LR1 for more Nanopore QC. Following the report, reads can be extracted using the `KrakenExtract` flag.   

Steps to run this pipeline:
1. Install a recent version of [nextflow](https://github.com/nextflow-io/nextflow).
2. Clone this repository using `git clone` into your working directory.
3. Create the conda environment using `conda env create --name ReadsQC --file ./Conda_Environments/ReadsQC.yml` from the Conda_Environment subsidectory [anaconda environment built: 11/2022]
4. Update the project directory path and add required (and optional) arguments in the `ReadsQC.sh` script, including the read type [`IlluminaSE`, `IlluminaPE`, or `Nanopore`].
5. Run the pipeline using `qsub ReadsQC.sh`.

For a more detailed help message run `sh ReadsQC.sh --help`

## Long Read Analyses

### LR1. Nanopore basecalling and fast5 QC
This pipeline extracts fast5 files from a tarball and uses [ONT-Guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) to call bases and emits a single fastq with all passed reads. Quality control plots are achieved through the use of [nanoQC](https://github.com/wdecoster/nanoQC) and [pycoQ](https://github.com/a-slide/pycoQC). If length filtering is needed, tools like seqtk can be used to filter by length, for example. 

Steps to run this pipeline:
1. Obtain the ONT-Guppy singularity container from me, or build the GPU version of [ONT-Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revaf_14dec2018/linux-guppy) into a singularity container called `ONT_Guppy_GPU.sif`.
2. Install a recent version of [nextflow](https://github.com/nextflow-io/nextflow).
3. Download a copy of the appropriate reference genome (i.e., the *Anopheles gambiae* reference genome from [vectorbase](https://vectorbase.org/vectorbase/app/record/dataset/DS_2251b21396)). 
4. Clone this repository using `git clone` into your working directory.
5. Create the conda environment using `conda env create --name NanoQC --file ./Conda_Environments/NanoQC.yml` from the Conda_Environment subsidectory [anaconda environment built: 11/2022]
6. Update the project directory path, ephemeral mount path, export the nextflow program path, and add required (and optional) arguments in the `ONT_Basecalling_QC.sh` script.
7. Run the pipeline using `qsub ONT_Basecalling_QC.sh`.

For a more detailed help message run `sh ONT_Basecalling_QC.sh --help`

### LR2. Consensus reference guided and *de novo* isoform transcriptome analysis 
This pipeline was developed to map full transcipts and use *de novo* transcript isoform assembly tools. The tools were tested using ONT cDNA libraries and are not appropriate for illumina data, but will work for PacBio isoseq and ONT dRNA approaches. This pipeline requires either filtered fastq files or bam files as input. If fastq files are provided, reads are mapped to a reference genome using either [uLTRA](https://github.com/ksahlin/ultra) or [MiniMap2](https://github.com/lh3/minimap2), with the former being considerably more computationally demanding. 

A consensus based approach is applied by default where novel isoforms are called using [IsoQuant](https://github.com/ablab/IsoQuant) and [ESPRESSO](https://github.com/Xinglab/espresso). A consensus call is created among tools using [GFFcompare](https://github.com/gpertea/gffcompare), but all gtfs will be emitted in case of high discordance among calls.

Steps to run this pipeline:
1. Install a recent version of [nextflow](https://github.com/nextflow-io/nextflow).
2. Clone this repository using `git clone` into your working directory.
3. Create the conda environment using `conda env create --name NanoIsoExpress --file ./Conda_Environments/NanoIsoExpress.yml` from the Conda_Environment subsidectory [anaconda environment built: 02/2023].
4. Clone the [ESPRESSO](https://github.com/Xinglab/espresso/tree/main/src) source perl scripts and change permissions to enable execution.  
5. Update the project directory path and add required (and optional) arguments in the `RNA_Isoform.sh` script.
6. Run the pipeline using `qsub RNA_Isoform.sh`.

For a more detailed help message run `sh RNA_Isoform.sh --help`. Also, please note that if fastq files are used as input, decompressed fastqs are needed for uLTRA. 


 
