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

  tag { id }
  label (params.LABEL)

  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(id), path(fasta)
  path(filterfasta)

  output:
  tuple val(id), path("${id}.cln"), emit: cln
  tuple val(id), path("${id}.cidx"), emit: cidx
  tuple val(id), path("${id}.clean"), emit: clean

  """
  export USER=\$(id -u -n)
  /usr/local/src/PASApipeline/bin/seqclean ${fasta} -v ${filterfasta} -c ${task.cpus} ${params.EXTRAPARS}
  """

}

// Consider including droping database via parameter
process importMySQLPasa {

  label (params.LABEL)

  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

	input:
	path(pasaconffilegeneral)
  path(dbparams)
  path(pasaschema)

	output:
	path("conftxt.new")

	"""
  source ${dbparams}
	mysql -u\${dbuser} -p\${dbpass} -h\${dbhost} -P\${dbport} -e "DROP DATABASE IF EXISTS \${dbname}; CREATE DATABASE \${dbname};"
	mysql -u\${dbuser} -p\${dbpass} \${dbname} -h\${dbhost} -P\${dbport} < ${pasaschema} > done_mysql

 	# Simple modification. This would need more love
 	cp ${pasaconffilegeneral} conftxt.new
	sed -i '/^MYSQLSERVER=/d' conftxt.new
	echo "MYSQLSERVER=\${dbhost}:\${dbport}" >> conftxt.new
	"""

}

process runPASA {

  tag { id }
  label (params.LABEL)

  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(id), path(relatedfasta)
  path(genome)
  path(seqclean_idx)
  path(seqclean_clean)
  path(seqclean_cln)
  path(model_gtf_file)
  path(pasaconffiledb)
  path(conftxt)
  val(pasamode)

  output:
  path("pasa*"), emit: pasa_files
  path("*assemblies.fasta"), emit: pasa_fasta
  path("*assemblies.gff3"), emit: pasa_gff3

	"""
 /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -PASACONF ${conftxt} -c ${pasaconffiledb} \
 -R -g ${genome} -t ${seqclean_clean} -T -u ${relatedfasta} \
 --trans_gtf ${model_gtf_file} --ALIGNERS ${pasamode} --TRANSDECODER --CPU ${task.cpus}
	"""

}


process generatePASAtrainingSet {

  tag { pasa_assemblies_fasta }
  label (params.LABEL)

  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  path(pasa_assemblies_fasta)
  path(pasa_assemblies_gff3)
  path(conftxt)

  output:
  path("*.transdecoder.gff3"), emit: transdecoder_gff3
  path("*.transdecoder.pep"), emit: transdecoder_pep

  """
  /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --PASACONF ${conftxt} --pasa_transcripts_fasta ${pasa_assemblies_fasta} --pasa_transcripts_gff3 ${pasa_assemblies_gff3} --single_best_only
  """

}


workflow PASA_SEQ_CLEAN {
    take:
    fasta
    filterfasta

    main:
    seqClean(fasta, filterfasta)

    emit:
    cln = seqClean.out.cln
    cidx =  seqClean.out.cidx
    clean =  seqClean.out.clean

}

workflow PASA_IMPORT_MYSQL {
    take:
    pasaconffilegeneral
    dbparams
    pasaschema

    main:
    out = importMySQLPasa(pasaconffilegeneral, dbparams, pasaschema)

    emit:
    out

}

workflow PASA_RUN_PASA {
    take:
    relatedfasta
    genome
    seqclean_idx
    seqclean_clean
    seqclean_cln
    model_gtf_file
    pasaconffiledb
    conftxt
    pasamode

    main:
    out = runPASA(relatedfasta, genome, seqclean_idx, seqclean_clean, seqclean_cln, model_gtf_file, pasaconffiledb, conftxt, pasamode)

    emit:
    out.pasa_files
    out.pasa_fasta
    out.pasa_gff3

}

workflow PASA_GENERATE_TRAINING_SET {
    take:
    pasa_assemblies_fasta
    pasa_assemblies_gff3
    conftxt

    main:
    out = generatePASAtrainingSet(pasa_assemblies_fasta, pasa_assemblies_gff3, conftxt)

    emit:
    out.transdecoder_gff3
    out.transdecoder_pep

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}
