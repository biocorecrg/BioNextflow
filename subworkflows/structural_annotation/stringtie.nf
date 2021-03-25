/*
* Stringtie subworkflow
* The parameters are:
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/stringtie:2.0--hc900ff6_0"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    echo stringtie' '`stringtie --version`
    """
}

process runStringtie {

  tag { id }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(id), path(bam_alignment)

  output:
  tuple val(id), path("${id}.stringtie.gtf")

  """
  stringtie ${bam_alignment} ${params.EXTRAPARS} --rf -p ${task.cpus} -v -o ${id}.stringtie.gtf
  """

}

workflow STRINGTIE {
    take:
    bam_alignment

    main:
    out = runStringtie(bam_alignment)

    emit:
    out

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}
