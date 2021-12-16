/*
* kraken2 subworkflow
* The parameters are:
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = "kraken2"
params.CONTAINER = ""

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    echo '2 latest version'
    """
}

process kraken2_build {

  tag { id }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:

  output:

  script:
  """
  """

}

process kraken2 {

  tag { id }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(pair_id), path(reads)

  output:
  path("kraken2*.report")
  path("kraken2*.out")
  path("cfs*.fq.gz")
  path("ucfs*.fq.gz")

  script:
  """
  """

}

workflow KRAKEN2_BUILD {
    take:


    main:

    out = kraken2_build()

    emit:
    out

}

workflow KRAKEN2 {
    take:


    main:

    out = kraken2()

    emit:
    out

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}
