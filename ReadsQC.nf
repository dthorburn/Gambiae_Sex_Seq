#!/usr/bin/env nextflow

/*
 *  Pipeline developed for the QC of raw reads. This pipeline was specifically designed with
 *  illumina short read and nanopore long reads in mind.  
 *  Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
 *  Date last modified: 30/11/2022
 */
                                                            // ========================================================
                                                            // Setting the help messages up
                                                            // ========================================================
def helpMessage() {
  log.info """
Rationale:
These pipelines were developed for a more thorough raw short and long read quality control measures.
They can handle Illumina paired-end and single-end data, and nanopore fastq files too.

ReadTypes:
[IlluminaPE/SE]:  FastP, FastQC, FastQ-Screen, and Kraken2 --> MultiQC
[Nanopore]:       FastP                        and Kraken2 --> MultiQC

Usage:
1. Create the conda3 environment using 'conda env create --name ReadsQC --file ReadsQC.yml' 
2. Update project directory in ReadsQC.sh
3. Add required arguments below
4. Submit pipeline coordinator using 'qsub ReadsQC.sh'

If required, change the concurrency options to free up jobs for other tasks. ICL HPC limits users to
50 jobs max. If you breach this limit, the coordinator may fail. 

Required arguments:
  --InDir                                     Path to directory with fastq files (ideally, gzip compressed files)
  --ReadType                                  ReadType selection [IlluminaSE, IlluminaPE, Nanopore, and KrakenExtract]
  --KrakenDB                                  Path to Kraken2 database

Concurrency arguments:
  --FP_Forks                                  Number of concurrent FastP jobs. Default: 10
  --KN_Forks                                  Number of concurrent Kraken2 jobs. Default: 10
  --FQC_Forks                                 Number of concurrent FastQC jobs. Default: 10
  --FQS_Forks                                 Number of concurrent Fastq-Screen jobs. Default: 10
  --MQS_Forks                                 Number of concurrent MultiQC jobs. Default: 10
  --KNE_Forks                                 Number of concurrent Kraken-Extract jobs. Default: 10

Additional arguments:
  --KrakenExtract                             Extract reads tagged with NCBI ID from raw reads (No other processes run)
  --KNE_TaxID                                 NCBI taxonomic ID used to extract reads from raw reads (i.e., 7165 for A. gambiae; 7164 
                                              for the genus Anopheles)
  --FP_threads                                Number of threads for each subprocess - swap FP for other jobs as needed (i.e., KN, FQS,
                                              FQC, MQC, or KNE)
  --FP_memory                                 Number of Gb of memory for each subprocess 
  --FP_walltime                               Number of hours for each subprocess (72 is maximum)
  --FP_args                                   Additional arguments for FastP - again, swap FP for other jobs as needed
  --FQS_conf                                  Path to template FastQ-Screen configuration file (includes path to FQS databases)
  --FQS_aligner                               Which aligner used by FastQ-Screen (bowtie, bowtie2, or bwa; default: bowtie2)
======================================================================================================================================
  """
}

log.info """
======================================================================================================================================
                                              ONT and Illumina QC Pipeline v0.1
======================================================================================================================================
Project Directory           : ${params.ProjectDIR}
KrakenDB                    : ${params.KrakenDB}
Raw Reads                   : ${params.InDir}
ReadType                    : ${params.ReadType}
Raw QC Reports              : ${params.RawQCDir}
MultiQC                     : ${params.MQCDir}
Kraken2 Filtered Reads      : ${params.KNEDir}
Filtering Reads by TaxID    : ${params.KrakenExtract}

Pipeline Author             : DMJ Thorburn <d.thorburn@imperial.ac.uk> - Please cite dthorburn/Gambiae_Sex_Seq GitHub repo in methods!
======================================================================================================================================
"""

if (params.help) {
    helpMessage()
    exit 0
}

if(!params.InDir) {
  log.info"""
ERROR: No input path to raw reads directory provided! --Input /path/to/fastqs/
======================================================================================================================================
  """
  helpMessage()
  exit 0
}
                                                            // ========================================================
                                                            // Setting up the channels
                                                            // ========================================================
// =================
// Value channels 
// =================
Channel
  .fromPath("${params.KrakenDB}")
  .ifEmpty { error "ERROR: No Kraken2 database path provided! --KrakenDB /path/to/krakenDB/" }
  .set { KDB }

Channel
  .fromPath("${params.ProjectDIR}")
  .ifEmpty { error "ERROR: No Kraken2 database path provided! --KrakenDB /path/to/krakenDB/" }
  .set { PDIR }

KrakenDB         = KDB.first()                                    // Unsure if there is a more appropriate way of setting up a value channel that is a path. 
projectDIR       = PDIR.first()   
ReadType         = Channel.value( params.ReadType )
conf_template    = file( params.FQS_conf, checkIfExists: true )


// =================
// Queue channels 
// =================
// Setting up the multiple parallel input channels for QC
if( params.ReadType == "IlluminaSE" || params.ReadType == "Nanopore" ){
  Channel
    .fromPath("${params.InDir}/*{.fastq.gz,.fq.gz,.fastq,.fq}")
    .ifEmpty { error "Cannot find any SE fastq files in ${params.InDir}" }
    .map { it -> tuple( it.simpleName, it ) }
    .into { reads_fastp; reads_fastqscreen; reads_fastqc; reads_kraken; reads_kne }
} else if( params.ReadType == "IlluminaPE" ){
  Channel
    .fromFilePairs("${params.InDir}/*_{R1,R2,1,2}{.fastq.gz,.fq.gz,.fastq,.fq,_001.fastq.gz,_001.fq.gz,_001.fastq,_001.fq}")
    .ifEmpty { error "Cannot find any PE fastq files in ${params.InDir}" }
    .into { reads_fastp; reads_fastqscreen; reads_fastqc; reads_kraken; reads_kne }
//} else if( params.ReadType == "Nanopore" ){
  // Even though they are not used, a channel for FQC and FQS are still needed.  
//  Channel
//    .fromPath("${params.InDir}/*{.fastq.gz,.fq.gz,.fastq,.fq}")
//    .ifEmpty { error "Cannot find any fastq files in ${params.InDir}" }
//    .into { reads_fastp; reads_fastqscreen; reads_fastqc; reads_kraken; reads_kne }
} else {
  log.info"""
ERROR: ReadType "${params.ReadType}" not recognised. Use either 'IlluminaPE', 'IlluminaSE', 'Nanopore', or 'KrakenExtract'
"""
helpMessage()
exit 0
}

// Separated the ReadType (previously Mode) and KrakenExtract as I needed the read type to handle the number of read files correctly in KNE.  
if( params.KrakenExtract ){
  Channel
    .fromPath("${params.RawQCDir}/kraken/*.out.txt")
    .ifEmpty { error "Cannot find the kraken output and report files in ${params.RawQCDir}/kraken/" }
    .map { file -> tuple(file.baseName.replaceAll(".out",""), file) }
    .set{ report_kne }

  reads_kne
    .combine( report_kne, by : 0 )
    .set { kne_in }
} else {
  // I don't understand why this is needed when the conditions of the process are not met but I get an error if I don't have a kne_in channel.
  // Mimicking the structure and type of input block too as that's also important... 
  Channel
    .fromPath( "${params.InDir}/*{.fastq.gz,.fq.gz,.fastq,.fq}" )
    .first()
    .map { it -> tuple(it.baseName, it, it.getParent()) }
    .set{ kne_in }
}
                                                            // ======================
                                                            // Step 1: FastP
                                                            // ======================

process FastP {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  maxForks params.FP_Forks
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate ReadsQC'
  
  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.FP_threads}:mem=${params.FP_memory}gb -lwalltime=${params.FP_walltime}:00:00"
  tag{ sampleID }

  publishDir(
    path: "${params.RawQCDir}/fastp",
    mode: 'copy',
  )

  input:
  val ReadType
  set sampleID, path(reads) from reads_fastp


  output:
  tuple sampleID, path("*.json"), path("*.html") into FP_Out

  when:
  ( ReadType =~ /^[IN]/ && !params.KrakenExtract )

  // Handling both cases: 
  // if(ReadType == "IlluminaPE") tuple == [sampleID, [readsR1.fq, readsR2.fq]]
  // if(ReadType != "IlluminaPE") tuple == [sampleID, reads.fq]
  script:
  if( params.ReadType == "IlluminaPE" ){
    """
    n_slots=`expr ${params.FP_threads} - 1`
    taskset -c 0-\${n_slots} fastp -i ${reads[0]} -I ${reads[1]} --json ${sampleID}_fastp.json --html ${sampleID}_fastp.html --thread ${params.FP_threads} ${params.FP_args}
    """
  } else {
    """
    n_slots=`expr ${params.FP_threads} - 1`
    taskset -c 0-\${n_slots} fastp -i ${reads} --json ${sampleID}_fastp.json --html ${sampleID}_fastp.html --thread ${params.FP_threads} ${params.FP_args}
    """
  }
}

                                                            // ======================
                                                            // Step 2: Kraken2
                                                            // ======================

process Kraken2 {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  maxForks params.KN_Forks
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate ReadsQC'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.KN_threads}:mem=${params.KN_memory}gb -lwalltime=${params.KN_walltime}:00:00"
  tag{ sampleID }

  publishDir(
    path: "${params.RawQCDir}/kraken",
    mode: 'copy',
  )

  input:
  val ReadType
  path KrakenDB
  set sampleID, path(reads) from reads_kraken

  output:
  tuple sampleID, path("*.report.txt"), path("*.out.txt") into KN_Out

  when:
  (ReadType =~ /^[IN]/ && !params.KrakenExtract )

  script:
  if( params.ReadType == "IlluminaPE" ){
    """
    ## Handling when files are not gzipped. If files are bgzipped, this will likely fail. Clunky, but it works. 
    if [[ ${reads[0]} =~ ".gz" ]]
    then
      echo "Input files compressed"
      cp ${reads[0]} R1.fq.gz
      cp ${reads[1]} R2.fq.gz
    else
      echo "Compressing files"
      gzip ${reads[0]} > R1.fq.gz
      gzip ${reads[1]} > R2.fq.gz
    fi

    kraken2 \\cd ../../
      --db ${KrakenDB} \\
      --threads ${params.KN_threads} \\
      --paired \\
      --report ${sampleID}.report.txt \\
      --output ${sampleID}.out.txt \\
      --gzip-compressed ${params.KN_args} \\
      R1.fq.gz \\
      R2.fq.gz
    """
  } else {
    """
    if [[ ${reads} =~ ".gz" ]]
    then
      echo "Input files compressed"
      cp ${reads} Reads.fq.gz
    else
      echo "Compressing files"
      cp ${reads} Reads.fq.gz
    fi

    kraken2 \\
      --db ${KrakenDB} \\
      --threads ${params.KN_threads} \\
      --report ${sampleID}.report.txt \\
      --output ${sampleID}.out.txt \\
      --gzip-compressed \\
      Reads.fq.gz ${params.KN_args}
    """ 
  }
}

                                                            // ======================
                                                            // Step 3: FastQC
                                                            // ======================
process FastQC {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  maxForks params.FQC_Forks
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate ReadsQC'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.FQC_threads}:mem=${params.FQC_memory}gb -lwalltime=${params.FQC_walltime}:00:00"
  tag{ sampleID }

  publishDir(
    path: "${params.RawQCDir}/fastqc",
    mode: 'copy',
  )

  input:
  val ReadType
  set sampleID, path(reads) from reads_fastqc

  output:
  tuple sampleID, path("*.html"), path ("*.zip") into FQC_Out

  when:
  (ReadType =~ /llumina/ && !params.KrakenExtract )

  script:
  // Handling when it's both SE and PE ReadTypes for Illumina reads
  // Will collect break this and do everything in every job?
  def readsQC = reads.collect{ "fastqc -t ${params.FQC_threads} -d ./tmp/ ${params.FQC_args} $it" }.join('\n')
  """
  mkdir tmp
  ${readsQC}
  """
}

                                                            // ======================
                                                            // Step 4: FastQ-Screen
                                                            // ======================
process FastQ_Screen {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  maxForks params.FQS_Forks
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate ReadsQC'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.FQS_threads}:mem=${params.FQS_memory}gb -lwalltime=${params.FQS_walltime}:00:00"
  tag{ sampleID }

  publishDir(
    path: "${params.RawQCDir}/fastqscreen",
    mode: 'copy',
  )

  input:
  val ReadType
  path conf_template
  set sampleID, path(reads) from reads_fastqscreen


  output:
  tuple sampleID, path("*.html"), path("*.png"), path("*.txt") into FQS_Out

  when:
  (ReadType =~ /llumina/ && !params.KrakenExtract )

  script:
  """
  ## Editting the configuration files
  cp ${conf_template} fastq_screen.conf
  bwt_path=`which bowtie`
  bwt2_path=`which bowtie2`
  bwa_path=`which bwa`
  sed -i "s|bwt_edit|\${bwt_path}|" fastq_screen.conf
  sed -i "s|bwt2_edit|\${bwt2_path}|" fastq_screen.conf
  sed -i "s|bwa_edit|\${bwa_path}|" fastq_screen.conf
  sed -i "s|edit_me|${params.FQS_threads}|" fastq_screen.conf

  fastq_screen --conf fastq_screen.conf --aligner ${params.FQS_aligner} ${params.FQS_args} ${reads}
  """
}
                                                            // ======================
                                                            // Step 5: MultiQC
                                                            // ======================
// A method to group outputs and regardless of the number of input channels, the resultant queue channel 
// should have the same file structure permitting a common input declaration for MQC. 
// [sampleID, [FQS_Files, KN_Files, FP_Files, FQC_Files]] - IlluminaPE/SE  
// [sampleID, [KN_Files, FP_Files]]                       - Nanopore


// Ended up going with a simpler solution. Once all processes are collected and complete, a channel will be generated with only the first
// sampleID omitted. The process will instead use the path to the QC directory as its input. A lazy work around, but it works. 
if( ReadType =~ /llumina/ && !params.KrakenExtract ){
  FQS_Out
    .mix(FQC_Out, KN_Out, FP_Out)
    .collect()
    .flatten()
    .first()
    .set{ MQC_In }
} else if( ReadType =~ /anopore/ && !params.KrakenExtract ){
  KN_Out
    .join(FP_Out)
    .collect()
    .flatten()
    .first()
    .set{ MQC_In }
} else {
  // I don't understand why this is needed when the conditions of the process are not met but I get an error if I don't have a MQC_In channel.  
  Channel
    .fromList( ['whymustyoutormentmenextflow'] )
    .set{ MQC_In }
}

process MultiQC {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  maxForks params.MQC_Forks
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate ReadsQC'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.MQC_threads}:mem=${params.MQC_memory}gb -lwalltime=${params.MQC_walltime}:00:00"

  publishDir(
    path: "${params.MQCDir}",
    mode: 'copy',
  )

  input:
  val ReadType
  path projectDIR
  val sampleID from MQC_In

  when:
  ( ReadType =~ /^[IN]/ && !params.KrakenExtract )

  output:
  path("*.html")


  script:
  if( params.ReadType == "Nanopore" ){
    """
    multiqc ${params.MQC_args} \\
      ${projectDIR}/${params.RawQCDir}/fastp \\
      ${projectDIR}/${params.RawQCDir}/kraken 
    """
  } else {
    """
    multiqc ${params.MQC_args} \\
      ${projectDIR}/${params.RawQCDir}/fastp \\
      ${projectDIR}/${params.RawQCDir}/fastqc \\
      ${projectDIR}/${params.RawQCDir}/fastqscreen \\
      ${projectDIR}/${params.RawQCDir}/kraken
    """
  }
}
                                                            // ======================
                                                            // Step 6: Kraken Extract
                                                            // ======================
process KrakenExtract {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  maxForks params.KNE_Forks
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate ReadsQC'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.KNE_threads}:mem=${params.FQS_memory}gb -lwalltime=${params.FQS_walltime}:00:00"

  publishDir(
    path: "${params.KNEDir}",
    mode: 'copy',
  )

  input:
  val ReadType
  set sampleID, path(reads), path(KN_out) from kne_in

  output:
  path("*.fastq.gz")

  when:
  params.KrakenExtract

  script:
  if( params.ReadType == "IlluminaPE" ){
    """
    py_script=`which extract_kraken_reads.py`
    \${py_script} \\
      -k ${KN_out} \\
      -s ${reads[0]} \\
      -s2 ${reads[1]} \\
      -o ${sampleID}_NCBI${params.KNE_TaxID}_R1.fastq \\
      -o2 ${sampleID}_NCBI${params.KNE_TaxID}_R2.fastq \\
      -t ${params.KNE_TaxID} \\
      --fastq-output ${params.KNE_args}

    gzip ${sampleID}_NCBI${params.KNE_TaxID}_R1.fastq
    gzip ${sampleID}_NCBI${params.KNE_TaxID}_R2.fastq
    """
  } else {
    """
    py_script=`which extract_kraken_reads.py`
    \${py_script} \\
      -k ${KN_out} \\
      -s ${reads} \\
      -o ${sampleID}_NCBI${params.KNE_TaxID}.fastq \\
      -t ${params.KNE_TaxID} \\
      --fastq-output ${params.KNE_args}

    gzip ${sampleID}_NCBI${params.KNE_TaxID}.fastq
    """
  }
}
