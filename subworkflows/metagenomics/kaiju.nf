/*
* kaiju subworkflow
* The parameters are:
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/kaiju:1.9.2--h5b5514e_0"

include { separateSEandPE } from '../global_functions.nf'


process getVersion {

    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    kaiju -h 2>&1| grep Kai 
    """
}

process alnSE {

  tag { pair_id }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(pair_id), path(reads)
  path(index)

  output:
  tuple val(pair_id), path("${pair_id}.kaiju"), emit: output

  script:
  """
  kaijux -z ${task.cpus} ${params.EXTRAPARS} -f ${index} -v -i ${reads} -o ${pair_id}.kaiju
  """

}

process alnPE {

  tag { pair_id }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(pair_id), path(readsA), path(readsB)
  path(index)

  output:
  tuple val(pair_id), path("${pair_id}.kaiju"), emit: output

  script:
  """
  kaijux -z ${task.cpus} ${params.EXTRAPARS} -f ${index} -v -i ${readsA} -j ${readsB} -o ${pair_id}.kaiju
  """

}

workflow ALN {

    take:
    fastq
    index

    main:
    def sep_fastq = separateSEandPE(fastq)
    sep_fastq.pe.view()
    outpe = alnSE(sep_fastq.se, index).output
    outse = alnPE(sep_fastq.pe, index).output

    //emit:
  //  out = outpe.mix(outse)

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    getVersion.out
}
