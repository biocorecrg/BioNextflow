/*
* Strawberry subworkflow
* The parameters are:
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "ruolinliu/strawberry:v1.1.1"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    echo `/home/strawberry/bin/strawberry 2>&1 |head -n2|tail -n1`
    """
}

process runStrawberry {

  tag { strawberry }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(id), path(bam_alignment)

  output:
  tuple val(id), path("${id}.strawberry.gtf")

  """
  /home/strawberry/bin/strawberry -p ${task.cpus} ${params.EXTRAPARS} -o ${id}.strawberry.gtf ${bam_alignment}
  """

}

workflow STRAWBERRY {
    take:
    bam_alignment

    main:
    out = runStrawberry(bam_alignment)

    emit:
    out

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}
