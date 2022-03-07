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
params.CONTAINER = "biocorecrg/kraken2:202112"

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

  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  val(groups)
  val(dbname)

  output:
  path(dbname)

  script:
  """
  kraken2-build --download-taxonomy --db ${dbname}
  orgs=${groups}
  for o in \${orgs//,/ }
  do
          kraken2-build --download-library \$o --db ${dbname}
          sleep 30
  done
  kraken2-build --build --db ${dbname}

  """

}

process kraken2 {

  tag { pair_id }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(pair_id), path(reads)
  path(database)

  output:
  path("kraken2*.report"), emit: report
  path("kraken2*.out"), emit: output
  path("cfs*.fq.gz"), emit: classified
  path("ucfs*.fq.gz"), emit: unclassified

  script:
  """
  mode=""
  if [[ "${reads}" = *" "* ]]; then
    mode="--paired"
  fi
  kraken2 --db ${database} --report kraken2_${pair_id}.report --threads ${task.cpus} \${mode} ${reads} --classified-out cfs_${pair_id}#.fq --unclassified-out ucfs_${pair_id}#.fq ${params.EXTRAPARS} > kraken2_${pair_id}.out
  gzip *.fq
  """

}

workflow BUILD {
    take:
    groups
    dbname

    main:

    out = kraken2_build(groups, dbname)

    emit:
    out

}


workflow RUN {
    take:
    fastq
    database

    main:
    out = kraken2(fastq, database)

    emit:
    out.report
    out.output
    out.classified
    out.unclassified

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    getVersion.out
}
