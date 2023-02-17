#!/usr/bin/env nextflow

/*
 *  Pipeline developed for calling bases using Oxford Nanopore fast5 files
 *  and for the QC of the subsequent calls. 
 *  Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
 *  Date last modified: 26/07/2022
 */
                                                            // ========================================================
                                                            // Setting the help messages up
                                                            // ========================================================
def helpMessage() {
  log.info """
  """
}

log.info """
==============================================================================================================================
                                              ONT-Guppy Basecalling and QC Pipeline v1
==============================================================================================================================
Guppy Container       : ${params.Container}
Guppy Config File     : ${params.Guppy_gupConf}
Tarball Directory     : ${params.Input}
Fastqs                : ${PWD}/02_Basecalled
QC Reports            : ${PWD}/03_QC
==============================================================================================================================
"""

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

if(!params.Container) {
  log.info"""
ERROR: No singularity container path provided! --Container /path/to/container.sif
==============================================================================================================================
  """
  helpMessage()
  exit 0
}

if(!params.Input) {
  log.info"""
ERROR: No input path to tarball directory provided! --Input /path/to/tarballs/
==============================================================================================================================
  """
  helpMessage()
  exit 0
}
                                                            // ========================================================
                                                            // Setting the value channels (can be read unlimited times)
                                                            // ========================================================
Guppy_container = file( params.Container, checkIfExists: true )
                                                            // ========================================================
                                                            // Step 1: Decompress tarball
                                                            // ========================================================
if ( params.Skip_Decompress == false ) {
  Channel
    .fromPath("${params.Input}/*.tar.gz")
    .ifEmpty { error "No tar.gz tarballs found in ${params.Input}!" }
    .map { file -> tuple(file.getSimpleName(), file) }
    .set { tarballs }
  // Had to do this for tb in ./*.tar.gz; do run_id=`echo "${tb}" | cut -d'_' -f 4`; echo $run_id; mv $tb ${run_id}.tar.gz; done
  process Decompress {
    errorStrategy { sleep(1200); return 'retry' }
    maxRetries 3
    cache = 'lenient'

    tag { run_id }

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.DC_threads}:mem=${params.DC_memory}gb -lwalltime=${params.DC_walltime}:00:00"

    publishDir(
      path: "${params.TarDir}",
      mode: 'move',
    )

    input:
    set val(run_id), path(tb) from tarballs

    output:
    tuple val(run_id), path("${run_id}/*.fast5") into fast5s
    path("${run_id}_inital_summary.txt")

    script:
    """
    tar -xf ${tb}
    mv ./*/final_summary.txt ${run_id}_inital_summary.txt
    mkdir ${run_id}
    mv ./*/fast5/*.fast5 ./${run_id}/
    """
  }
}

                                                            // ========================================================
                                                            // Step 2: ONT-Guppy Call
                                                            // ========================================================
if (params.Skip_Guppy == false) {
  if ( params.Skip_Decompress ){
    // Deprecated - cannot skip decompress at the moment. 
    //Channel
    //  .fromPath("${params.TarDir}/*.fast5")
    //  .ifEmpty { error "Cannot find fast5s in ${params.TarDir}" }
    //  .collect()
    //  .set { fast5s }
  }
  //Channel
  //  .fromPath("${params.TarDir}")
  //  .set { fast5_dir }
  //fast5s
  //  .map { run_id, files -> tuple(run_id, files, files.getParent()) }
  //  .set { fast5s_ch }

  process Guppy_Call { 
    errorStrategy { 'retry' }
    maxRetries 3
    cache = 'lenient'

    tag { run_id }
    
    publishDir(
      path: "${params.GuppyDir}",
      mode: 'copy',
    )
    
    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.Guppy_threads}:mem=${params.Guppy_memory}gb:ngpus=${params.Guppy_GPUs}:gpu_type=${params.Guppy_GPU} -lwalltime=${params.Guppy_walltime}:00:00"
   
    beforeScript 'module load cuda/11.4.2'
    
    input:
    set run_id, path(fast5_files) from fast5s
    path Guppy_container
    //path file_in from fast5_dir

    output:
    tuple val(run_id), path("${run_id}.txt") into summary_ch
    path("${run_id}.fastq.gz") into all_fqs

    script:
    """
    ##if [ ! -d indir ]; then mkdir indir; mv *.fast5 indir; fi
    if [ ! -d outdir ]; then mkdir outdir; fi
    export SINGULARITY_BIND="/rds/general/ephemeral/user/dthorbur/ephemeral/"
    n_slots=`expr ${params.Guppy_threads} / 2`
    taskset -c 0-\${n_slots} \\
      singularity exec --nv ${Guppy_container} \\
        guppy_basecaller \\
          --config "${params.Guppy_gupConf}" \\
          --device cuda:all \\
          --compress_fastq \\
          -i ${PWD}/01_Raw_Input/${run_id} \\
          -s outdir \\
          --gpu_runners_per_device "${params.Guppy_GPU_runners}" \\
          --num_callers "${params.Guppy_callers}" "${params.Guppy_args}"

    mv outdir/sequencing_summary.txt ${run_id}.txt
    cat outdir/pass/*.fastq.gz > ${run_id}.fastq.gz
    """
  }
}
                                                            // ========================================================
                                                            // Step 3: QC
                                                            // ========================================================
if (params.Skip_QC == false) {
  if (params.Skip_Guppy) {
    //Channel
    //  .fromPath("${params.GuppyDir}/*.fastq.gz")
    //  .ifEmpty { error "No fastq.gz found in ${params.GuppyDir}" }
    //  .set { all_fqs }
    Channel
      .fromPath("${params.GuppyDir}/*.txt")
      .ifEmpty { error "No summary txt reports found in ${params.GuppyDir}" }
      .set { summary_ch }
  }

  process NanoQC {
    errorStrategy { sleep(1200); return 'retry' }
    maxRetries 3
    cache = 'lenient'

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.QC_threads}:mem=${params.QC_memory} -lwalltime=${params.QC_walltime}:00:00"

    tag { run_id }

    publishDir(
      path: "${params.QCDir}",
      mode: 'copy',
    )

    beforeScript 'module load anaconda3/personal;source activate NanoPlot'

    input:
    set run_id, path(summary_file) from summary_ch

    output:
    path("${run_id}_pycoQC.html")
    path("${run_id}_NanoPlot.html")
    path("${run_id}_NanoPlot.tar.gz")

    script:
    """
    pycoQC --summary_file ${summary_file} --html_outfile ${run_id}_pycoQC.html && echo "Done with pycoQC: `date`"
    NanoPlot --summary ${summary_file} --loglength -o ${run_id}_NanoPlot && echo "Done with NanoPlot: `date`"

    cp ${run_id}_NanoPlot/NanoPlot-report.html ./${run_id}_NanoPlot.html
    tar -czvf ${run_id}_NanoPlot.tar.gz ${run_id}_NanoPlot/ && echo "Done compressing NanoPlot output: `date`"
    """ 
  }
}