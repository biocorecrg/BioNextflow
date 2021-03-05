/*
* CLASS subworkflow
* The parameters are:
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = ""

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    echo CLASS 0.0
    """
}

process runClass {

  tag { class }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(id), path(bam_alignment)

  output:
  tuple val(id), path("${id}.class.gtf")

  """
  perl /usr/local/src/CLASS/run_class.pl -a ${bam_alignment} -o ${id}.class.gtf -p ${task.cpus} --verbose
  """

}

workflow CLASS {
    take:
    bam_alignment

    main:
    out = runClass(bam_alignment)

    emit:
    out

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}
