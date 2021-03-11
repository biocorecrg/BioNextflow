/*
* Pasapipeline subworkflow
* The parameters are:
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "biocorecrg/pasapipeline:2.3.3"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    echo `/usr/local/src/PASApipeline/Launch_PASA_pipeline.pl --version`
    """
}


process seqClean {

  tag { pasapipeline }
  label (params.LABEL)

  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(id), path(fasta)
  path(filterfasta)

  output:
  tuple val(id), path("${id}.cln")
  tuple val(id), path("${id}.cidx")
  tuple val(id), path("${id}.clean")

  """
  export USER=\$(id -u -n)
  /usr/local/src/PASApipeline/bin/seqclean ${fasta} -v ${filterfasta} -c ${task.cpus} ${params.EXTRAPARS}
  """

}

// Consider including droping database via parameter
process importMySQLPasa {

  tag { pasapipeline }
  label (params.LABEL)

  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

	input:
	path(pasaconffilegeneral)

	output:
	path("done_mysql")
	path("conftxt")

	"""
	mysql -u$params.dbuser -p$params.dbpass -h$dbhost -P$params.dbport -e "DROP DATABASE IF EXISTS $params.dbname; CREATE DATABASE $params.dbname;"
	mysql -u$params.dbuser -p$params.dbpass $params.dbname -h$dbhost -P$params.dbport < $params.pasaschema > done_mysql

 	# Simple modification. This would need more love
 	cp ${pasaconffilegeneral} conftxt
	sed -i '/^MYSQLSERVER=/d' conftxt
	echo "MYSQLSERVER=${dbhost}:${params.dbport}" >> conftxt
	"""

}

process runPASA {

	label 'pasa'

 publishDir outputdir, mode: 'copy'

 input:
 file (relatedfasta) from file( params.relatedfasta )
 file (genome) from file( params.genome )
 file (relatedfasta_cidx) from seqclean_idx
 file (relatedfasta_clean) from seqclean_clean
 file (relatedfasta_cln) from seqclean_cln
 file (done_mysql) from done_mysql
 file (model_gtf_file) from model_gtf_file


 output:
 file "pasa*" into pasa_results
 file "*assemblies.fasta" into pasa_assemblies_fasta
 file "*assemblies.gff3" into pasa_assemblies_gff3


	// TODO: Docker and Singularity options. Bad for cloud. It seems it can be changed with --PASACONF

	containerOptions "--bind ${outputdir}/conftxt:/usr/local/src/PASApipeline/pasa_conf/conf.txt"

	"""
 /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c $params.pasaconffiledb \
 -R -g ${genome} -t ${relatedfasta_clean} -T -u ${relatedfasta} \
 --trans_gtf ${model_gtf_file} --ALIGNERS $params.pasamode --TRANSDECODER --CPU ${task.cpus}
	"""

}


process generatePASAtrainingSet {

 label 'pasa'

 publishDir outputdir+"/training", mode: 'copy'

 input:
 file (pasa_assemblies_fasta) from pasa_assemblies_fasta
 file (pasa_assemblies_gff3) from pasa_assemblies_gff3

 output:
 file "*.transdecoder.gff3" into pasa_transdecoder_training_gff3
 file "*.transdecoder.pep" into pasa_transdecoder_training_pep

 	// TODO: Docker and Singularity options. Bad for cloud
	containerOptions "--bind ${outputdir}/conftxt:/usr/local/src/PASApipeline/pasa_conf/conf.txt"

	"""
 /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta ${pasa_assemblies_fasta} --pasa_transcripts_gff3 ${pasa_assemblies_gff3} --single_best_only
 """

}


workflow PASAPIPELINE {
    take:
    fasta
    filterfasta

    main:
    out = seqClean(fasta, filterfasta)

    emit:
    out

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}
