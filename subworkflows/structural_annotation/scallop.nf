/*
* Scallop subworkflow
* The parameters are:
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/scallop:0.10.5--hf5e1fbb_0"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    echo scallop' '`scallop --version`
    """
}

process runScallop {

  tag { id }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(id), path(bam_alignment)

  output:
  tuple val(id), path("${id}.scallop.gtf")

  """
  scallop -i ${bam_alignment} ${params.EXTRAPARS} -o ${id}.scallop.gtf --verbose --library_type first
  """

}

workflow SCALLOP {
    take:
    bam_alignment

    main:
    out = runScallop(bam_alignment)

    emit:
    out

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}
