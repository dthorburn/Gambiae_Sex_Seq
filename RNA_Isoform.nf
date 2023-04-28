#!/usr/bin/env nextflow

/*
 * Pipeline developed for mapping bulk full transcript reads using uLTRA or MiniMap2,
 * followed by de novo isoform assembly and quantification using ESPRESSO and IsoQuant.   
 * Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
 * Date last modified: 07/03/2023
 */

// ========================
// Help and log messages
// ========================

def helpMessage() {
log.info """
Rationale:
This pipeline was developed to map full transcipts and use reference guided de novo transcript assembly tools. 
The tools were tested using ONT cDNA libraries and are not appropriate for illumina data. 

Setup:
1. Create the conda environment using 'conda -env create --name NanoIsoExpress --file NanoIsoExpress.yml'
2. Update project directory in RNA_Isoform.sh
3. Add required arguments to nextflow command in RNA_Isoform.sh
4. Submit pipeline coordinator using 'qsub RNA_Isoform.sh'

Usage:
The pipeline takes either quality filtered fastq files or sorted bam files:
nextflow run RNA_Isoform.nf -c RNA_Isoform.config --profile imperial --InDir /path/to/fastqs --Annots /path/to/gff
or
nextflow run RNA_Isoform.nf -c RNA_Isoform.config --profile imperial --MAPDir /path/to/bams --Annots /path/to/gtf.gz

Required arguments:
  --InDir                            Path to directory with fastq files [Not required if bam files are provided]
  --MAPDir                           Path to directory with bam files [Not required if fastq files are provided]. Bams need to
                                     be encoded in sam version 1.4 (i.e., CIGAR strings with +/= instead of M).
  --RefGen                           Path to reference genome uncompressed fasta file
  --Annots                           Path to annotations file [gff|gtf|gff.gz|gtf.gz]

Additional arguments:
   -w                                Path to [ephemeral|scratch] directory. If using uLTRA, several Tb of temp data will be
                                     generated
  --Mapper                           Read mapping tool selection [ultra | minimap2; default ultra]
  --MAP_Forks                        Concurrency option for mapping jobs [default: 20]
  --MAP_ultra_rt                     Read type for ultra [default: --ont] 
  --MAP_ultra_args                   Additional arguments for ultra
  --MAP_minimap2_ax                  Preset mapping minimap2 parameter [default: splice]
  --MAP_minimap2_args                Additional arguments for minimap2
  --ISQ_dt                           IsoQuant data type [default: nanopore]
  --ESC_xarg_args                    Additional arguments for xargs parallelising ESPRESSO_C.pl
  --ESQ_cutoff                       Number of perfectly mapped reads to call novel isoform [default: 4] 
  --Email                            If an email address is provided, an email will notify you of pipeline completion. Does not
                                     work if the pipeline coordinator job crashes.
  --clear_cache_on_success           If the pipeline runs successfully, it will delete the work directory [true|false; default: 
                                     true]
  
HPC exectutor arguments:             These options apply to all processes. Change MP to job to tune [ISQ, ESS, ESC, ESQ, or GFC]
  --MP_threads                       Number of cpus to request
  --MP_memory                        Number of GB of memory to request
  --MP_walltime                      Number of hours to request [max: 72] 
  --MP_args                          Additional arguments for each tool
================================================================================================================================
"""
}

// Ensuring the correct mapped reads directory is emitted in the log and output files 
if ( params.MAPDir == "02_Mapped" ) {
  params.MapPath = "${workflow.launchDir}/02_Mapped" 
} else {
  params.MapPath = "${params.MAPDir}"
}

start_date = new java.util.Date()   

log.info """
================================================================================================================================
                                        Long read RNA de novo isoform assembly pipeline v1
================================================================================================================================
Started       : ${start_date}
Mapper        : ${params.Mapper}
Reference     : ${params.RefGen}
Annotations   : ${params.Annots}
Fastqs        : ${params.InDir}
Mapped reads  : ${params.MapPath}
IsoQuant      : ${workflow.launchDir}/03_IsoQuant
ESPRESSO      : ${workflow.launchDir}/04_ESPRESSO
GFFCompare    : ${workflow.launchDir}/05_GFFCompare

Author: Doko-Miles Thorburn <d.thorbur@ic.ac.uk>
================================================================================================================================
"""

if (params.help) {
    helpMessage()
    exit 0
}


// ========================
// Required path checks 
// ========================

if(!params.RefGen) {
  log.info"""
ERROR: No reference genome path provided! --RefGen /path/to/genome.fasta
================================================================================================================================
  """
  helpMessage()
  exit 0
}

if(!params.Annots) {
  log.info"""
ERROR: No annotation file path provided! --RefGen /path/to/genome.[gff|gtf|gff.gz|gtf.gz]
================================================================================================================================
  """
  helpMessage()
  exit 0
}



// ========================
// Value channels 
// ========================

Channel
  .fromPath("${params.ESPPath}")
  .ifEmpty { error "ERROR: No path to ESPRESSO provided. --ESSPath /path/to/espresso/src" }
  .set { ESP }
// A few extra lines that permits case insensitivity in the input of the --Mapper argument
Channel
  .value( params.Mapper )
  .map { it -> it.toString().toLowerCase() }
  .set { MapperVal }
ESPath    = ESP.first() // There has got to be a more efficient way to do this. Also, unsure how NF will take using absolute paths
ref_genome = file( params.RefGen, checkIfExists: true )

// ========================
// Queue channels 
// ========================

// Assuming that if no --InDir is supplied, then it was bams that were used as initial input
if(!params.InDir){
  // Additionally creating the fastq_in channel as channels with the correct structure are needed even if the when argument
  // is not met.

/* Likely deprecated - will remove once testing completes
 *  Channel
 *    .fromFilePairs("${params.MAPDir}/*{bam,bai}") { file -> file.name.replaceAll(/.bam|.bai$/,'') }
 *    .ifEmpty { error "ERROR: Cannot find any indexed bam files in ${params.MAPDir}." }
 *    .map { ID, files -> tuple(files[0], files[1])}
 *    .set { bams_iq }
 * Channel
 *  .fromFilePairs("${params.MAPDir}/*{bam,bai}") { file -> file.name.replaceAll(/.bam|.bai$/,'') }
 *    .ifEmpty { error "ERROR: Cannot find any indexed bam files in ${params.MAPDir}." }
 *    .map { ID, files -> tuple(ID, files[0], files[1])}
 *    .set { bams_es }
 */

//  Channel
//    .fromFilePairs("${params.MAPDir}/*{bam,bai}") { file -> file.name.replaceAll(/.bam|.bai$/,'') }
//    .ifEmpty { error "ERROR: Cannot find any indexed bam files in ${params.MAPDir}." }
//    .map { ID, files -> tuple(ID, files[0], files[1])}
//    .set { bams_fmb }
  Channel
    .fromPath("${params.MAPDir}/*.bam")
    .ifEmpty { error "ERROR: Cannot find any bam files in ${params.MAPDir}." }
    .map { it -> tuple( it.simpleName, it ) }
    .set { bams_fmb }

  Channel
    .fromPath("${params.MAPDir}/*{.bam}")
    .map { it -> tuple( it.simpleName, it ) }
    .set { fastq_in }
} else {
  Channel
    .fromPath("${params.InDir}/*{.fastq.gz,.fq.gz,.fastq,.fq}")
    .ifEmpty { error "ERROR: Cannot find any fastq files in ${params.InDir}." }
    .map { it -> tuple( it.simpleName, it ) }
    .set { fastq_in }
}
Channel
  .fromPath("${params.Annots}")
  .ifEmpty { error "ERROR: No annotation file provided. --Annots /path/to/[gtf|gff](.gz)?" }
  .map { it -> tuple( it.simpleName, it ) }
  .set { in_gtf }


// ========================
// Step 0: Annotation Check 
// ========================

process GFF_Check {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 2
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate NanoIsoExpress'
  
  // Running locally within the submittion script job. Doesn't need to queue - it's small and quick.
  executor = 'local'

  input:
  set sampleID, path(input) from in_gtf

  output:
  // Not all tools need a gtf, but to make the most accurate comparison, all annotations will be in 
  // AGAT gtf format.
  path "${sampleID}_AGAT.gtf" into queue_c_gtf

  script:
  """
  if [[ ${input} =~ ".gz" ]]
  then
    echo "Decompressing file: `date`"
    gzip -d ${input} 
  fi

  if [[ `ls -1 ${sampleID}.*` =~ ".gff" ]]
  then
    echo "Converting gff to AGAT gtf: `date`"
    agat_convert_sp_gff2gtf.pl --gff ${sampleID}.gff --gtf ${sampleID}_AGAT.gtf
  elif [[ `ls -1 ${sampleID}.*` =~ ".gtf" ]]
  then
    rsync --links ${sampleID}.gtf ${sampleID}_AGAT.gtf
  else 
    echo "Check ${input} is correct format"
    echo "Accepted are .gtf, .gff, .gtf.gz, and .gff.gz"
    exit 1
  fi
  """
}
// It seems you cannot create a value channel as output of a process
// So this is converting the output queue channel to a value channel
queue_c_gtf
  .first()
  .set { ref_gtf }

// ========================
// Step 1: Mapping 
// ========================
if( params.InDir ) {
  process Mapping {
    errorStrategy { sleep(600); return 'retry' }
    maxRetries 2
    maxForks params.MAP_Forks
    cache 'lenient'
    beforeScript 'module load anaconda3/personal; source activate NanoIsoExpress'

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.MAP_threads}:mem=${params.MAP_memory}gb -lwalltime=${params.MAP_walltime}:00:00"
    tag{ sampleID+"-"+MapperVal }

    publishDir(
      path: "${params.MAPDir}",
      mode: 'copy',
    )

    input:
    path ref_genome
    path ref_gtf
    val MapperVal
    set sampleID, path(reads) from fastq_in

    output:
    //tuple sampleID, path("*_ES.bam"), path("*_ES.bam.bai") into bams_es
    //tuple sampleID, path("*_IQ.bam"), path("*_IQ.bam.bai") into bams_iq
    //path("*_ES.bam")     into bams_es
    //path("*_ES.bam.bai") into bais_es
    //path("*_IQ.bam")     into bams_iq
    //path("*_IQ.bam.bai") into bais_iq
    tuple sampleID, path("*.bam") into bams_fmb

    script:
    if( MapperVal == "ultra" )
      // It was just easier to index the reference genome in each sample rather than handling the number and diversity
      // of files created by uLTRA index. Issue compounded by problem below. 
      // For whatever reason, uLTRA can't seem to handle symlinks. Whilst Nextflow doesn't like using absolute paths, 
      // it seems I've had little choice but to use them here. 
      """
      echo "Indexing reference genome: `date`"
      mkdir references
      workDIR=`pwd`
      uLTRA index \\
        `readlink -f ${ref_genome}` \\
        `readlink -f ${ref_gtf}` \\
        \${workDIR}/references

      ## Annoyingly uLTRA can't seem to handle gzipped input fastqs... This is really slow and pointless...
      if [[ ${reads} =~ ".gz" ]]
      then
        echo "Decompressing input: `date`"
        gzip -cd ${reads} > input.fastq
      else 
        ## Since uLTRA can't handle symlinks it seemed easier to just copy the file -- only takes a few minutes
        rsync --copy-links --progress ${reads} input.fastq
      fi

      echo "Aligning reads: `date`"
      mkdir raw_output
      uLTRA align \\
        `readlink -f ${ref_genome}` \\
        \${workDIR}/input.fastq \\
        \${workDIR}/raw_output \\
        --t ${params.MAP_threads} \\
        --index \${workDIR}/references \\
        --prefix ${sampleID} ${params.MAP_ultra_rt} ${params.MAP_ultra_args}

      echo "Converting SAM files: `date`"
      samtools view -bh raw_output/${sampleID}.sam | samtools sort > ${sampleID}.bam
      """
    else if( MapperVal == "minimap2" )
      """
      minimap2 -ax ${params.MAP_minimap2_ax} --eqx ${params.MAP_minimap2_args} \\
        -t ${params.MAP_threads} \\
        ${ref_genome} \\
        ${reads} > ${sampleID}.sam
      echo "Converting SAM files: `date`"
      samtools view -bh ${sampleID}.sam | samtools sort > ${sampleID}.bam
      """
    else
      error "Invalid mapper provided: ${params.Mapper}, use either ultra or minimap2"
  }
}

// Added an additional step. When this was in the above mapping script it complicated the input when bams were supplied. 
process FormatBams {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate NanoIsoExpress'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.FMB_threads}:mem=${params.FMB_memory}gb -lwalltime=${params.FMB_walltime}:00:00"
  tag { sampleID }
  
  input:
  path ref_genome
  set sampleID, path(bam) from bams_fmb

  output: 
  path("${sampleID}_ES.bam")     into bams_es
  path("${sampleID}_ES.bam.bai") into bais_es
  path("${sampleID}_IQ.bam")     into bams_iq
  path("${sampleID}_IQ.bam.bai") into bais_iq

  script:
  """
  ##rsync --progress ${bam} ${sampleID}_IQ.bam
  ## Can symlinks work? - It looks like they can. Remove this comment if no probems arise. 
  rsync -l ${bam} ${sampleID}_IQ.bam

  ## Just in case. Don't think this is needed anymore.
  ## samtools faidx ${ref_genome}

  reformat.sh -eoom in=${sampleID}_IQ.bam out=${sampleID}_ES.bam sam=1.3 

  samtools index ${sampleID}_IQ.bam
  samtools index ${sampleID}_ES.bam
  """
}

// ========================
// Step 2: IsoQuant 
// ========================

process IsoQuant {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate NanoIsoExpress'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.ISQ_threads}:mem=${params.ISQ_memory}gb -lwalltime=${params.ISQ_walltime}:00:00"

  publishDir(
    path: "${params.ISQDir}",
    mode: 'copy',
  )

  input:
  path ref_genome
  path ref_gtf
  // The collect() function messed up tuple cardinality, but since none of these variables are needed in the script
  // we can get around this by just calling them all files
  // set sampleID, path(bams), path(bais) from bams_iq.collect
  // Separated 
  path(bams) from bams_iq.collect()
  path(bais) from bais_iq.collect()

  output:
  path "**.extended_annotation.gtf" into IQ_gtf
  tuple path("**.tsv"), path("**.bed"), path("**.transcript_models.gtf"), path("**.log")

  script:
  def list_bams = bams.collect{"ls -d -1 $it >> bam_list.txt" }.join('\n')
  """
  echo "Creating bam list: `date`"
  ##ls -d -1 $PWD/*.bam > bam_list.txt
  ${list_bams}
  
  echo "Starting IsoQuant: `date`"
  mkdir output
  isoquant.py --force \\
    --threads ${params.ISQ_threads} \\
    --reference ${ref_genome} \\
    --genedb ${ref_gtf} \\
    --bam_list bam_list.txt \\
    --read_group file_name \\
    --data_type ${params.ISQ_dt} \\
    -o ./output ${params.ISQ_args}
  """
}

// ========================
// Step 3: ESPRESSO 
// ========================

process ESPRESSO_S {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate NanoIsoExpress'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.ESS_threads}:mem=${params.ESS_memory}gb -lwalltime=${params.ESS_walltime}:00:00"

  input:
  path ref_genome
  path ref_gtf
  path ESPath
  //set sampleID, path(bams), path(bais) from bams_es.collect()
  // Couldn't get it working with the tuple as input cardinality was messing up.
  path(bams) from bams_es.collect()
  path(bais) from bais_es.collect()

  output:
  // Couldn't figure out what ESPRESSO_C needs individually to work, so I will output it all!
  path "working_dir" into ESS_out

  script:
  //def index_bams = bams.collect{ "samtools index $it" }.join('\n')
  def make_list = bams.collect{"bam_name=`basename $it | sed -e 's/.bam//'`; echo -e $it\\\t\$bam_name >> file_list.tsv" }.join('\n')
  """
  echo "Creating list: `date`"
  ${make_list}

  echo "Starting ESPRESSO_S: `date`"
  mkdir working_dir
  perl ${ESPath}/ESPRESSO_S.pl ${params.ESS_args} \\
    -T ${params.ESS_threads} \\
    -L file_list.tsv \\
    -F ${ref_genome} \\
    -A ${ref_gtf} \\
    -O ./working_dir
  """
}

process ESPRESSO_C {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate NanoIsoExpress'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.ESC_threads}:mem=${params.ESC_memory}gb -lwalltime=${params.ESC_walltime}:00:00"

  input:
  path ref_genome
  path ref_gtf
  path ESPath
  path working_dir from ESS_out

  output:
  path "working_dir" into ESC_out

  script:
  // Not the prettiest thing, and it will be horrifically memory hungry. 
  // But since I couldn't figure out how to parallelise this element since
  // It has to be in the same otuput directory and I don't know what global 
  // files are changed. The threads are rounded down and split equally among
  // the xargs jobs
  """
  num_folders=`ls -dq working_dir/[0-9]*/ | wc -l`
  num_threads=\$(( ${params.ESC_threads} / \${num_folders} ))

  printf %s\\\\n \$(seq 0 `expr \$num_folders - 1`) | xargs ${params.ESC_xarg_args} -n 1 -I {} \\
    perl ${ESPath}/ESPRESSO_C.pl ${params.ESC_args} \\
      -T \${num_threads} \\
      -I ${working_dir} \\
      -F ${ref_genome} \\
      -X {}
  """
}

process ESPRESSO_Q {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate NanoIsoExpress'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.ESQ_threads}:mem=${params.ESQ_memory}gb -lwalltime=${params.ESQ_walltime}:00:00"

  publishDir(
    path: "${params.ESSDir}",
    mode: 'copy',
  )

  input:
  path ref_genome
  path ref_gtf
  path ESPath
  path working_dir from ESC_out

  output:
  path "*.gtf" into ES_gtf
  tuple path("*.esp"), path("*.tsv"), path("*.list"), path("*.fa")

  script:
  """
  ## Unsure if the output directory is needed and could be simplified if no other outputs are generated.
  mkdir OUT_DIR
  perl ${ESPath}/ESPRESSO_Q.pl ${params.ESQ_args} \\
    -T ${params.ESQ_threads} \\
    -L ${working_dir}/*.updated \\
    -A ${ref_gtf} \\
    -V OUT_DIR/compatible_isoform.tsv \\
    --read_num_cutoff ${params.ESQ_cutoff}
  """
}

// ========================
// Step 4: GFF Compare 
// ========================

process GFFCompare {
  errorStrategy { sleep(600); return 'retry' }
  maxRetries 3
  cache 'lenient'
  beforeScript 'module load anaconda3/personal; source activate NanoIsoExpress'

  executor = 'pbspro'
  clusterOptions = "-lselect=1:ncpus=${params.GFC_threads}:mem=${params.GFC_memory}gb -lwalltime=${params.GFC_walltime}:00:00"

  publishDir(
    path: "${params.GFCDir}",
    mode: 'move',
  )

  input:
  path ref_gtf
  path IQG from IQ_gtf
  path EQG from ES_gtf

  output:
  path "*"

  script:
  """
  gffcompare ${params.GC_args} -r ${ref_gtf} ${IQG} -o "AgamP4_IsoQuant"
  gffcompare ${params.GC_args} -r ${ref_gtf} ${EQG} -o "AgamP4_ESPRESSO"
  gffcompare ${params.GC_args} -r ${EQG} ${IQG} -o "ESPRESSO_IsoQuant" -p "ES_IQ_TCONS"

  echo "Long read RNA isoform analysis completed on: `date` > RNA_Isoform_Versions.txt 
  conda --version >> RNA_Isoform_Versions.txt
  conda list >> RNA_Isoform_Versions.txt
  """
}

workflow.onComplete {
  log.info "Executing workflowcomplete"
  // Generate report for email if provided
  def subject = "ICL-RNA Isoform Pipeline Complete: ${workflow.runName}"
  def msg = """\
      Long read RNA de novo isoform assembly pipeline summary
      -------------------------------------------------------
      Nextflow version  : ${nextflow.version}
      Run Name          : ${workflow.runName}
      Started           : ${workflow.start}
      Completed         : ${workflow.complete}
      Duration          : ${workflow.duration}
      Success           : ${workflow.success}
      Project directory : ${workflow.launchDir}
      Work directory    : ${workflow.workDir}
      Exit Status       : ${workflow.exitStatus}
      """
      .stripIndent()
  
  if ( report.enabled == true ) {
    def attmnt = file("${workflow.launchDir}/${report.file}")
  } else {
    // Need to see what happens when no attachment is provided
    boolean attmnt = false
  }

  // Sends email to provided address on completion.  
  if ( params.Email ) {
    log.info "[ICL-RNA_Isoform] Trying to send email"
//    try {
      sendMail(
        to: ${params.Email}, 
        subject: 'RNA_Isoform.nf Pipeline Execution',
        attach: attmnt,
        body: msg)
      log.debug "[ICL-RNA_Isoform] Sent summary e-mail using sendmail"
//    } catch (all) {
//      log.debug "[ICL-RNA_Isoform] Sendmail failed."
//    }
  }

  // If the pipeline completed successfully and the clear_cache_on_success option is set to true [default] it will delete 
  // the work directory
  if ( workflow.success == true && params.clear_cache_on_success ) {
    log.info "Workflow ${workflow.runName} completed successfully. Removing ${workflow.workDir} directory."
    //file('./work/').deleteDir()
    file("${workflow.workDir}").deleteDir() // Is this better or quite a dangerous paramater to set? 
  }
}