/*
* geneidtrainer subworkflow
* The parameters are:
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = "geneidtrainer"
params.CONTAINER = "biocorecrg/geneidtrainer:1.4.5"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    /** TODO: Actual tool should release version **/
    """
    echo '1.4.5'
    """
}

process geneidtrainer {

  tag { id }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(species), path(gff), path(fasta), val(reduced), path(userdata)

  output:
  tuple val(species), path("${species}.geneid.optimized.param")

  script:
  """
  REDUCED=no
  USERDATA=""
  if [ ! -z "${reduced}" ]
  then
    REDUCED=${reduced}
  fi
  if [ ! -z "${userdata}" ]
  then
    USERDATA="--userdata ${userdata}"
  fi

  /scripts/geneidTRAINer4docker.pl -species ${species} -gff ${gff} -fastas ${fasta} -reduced \$REDUCED \$USERDATA -results ${params.OUTPUT} ${params.EXTRAPARS}
  """

}

workflow GENEIDTRAINER {
    take:
    species
    gff
    fasta
    reduced
    userdata

    main:

    out = geneidtrainer(species, gff, fasta, reduced, userdata)

    emit:
    out

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}
