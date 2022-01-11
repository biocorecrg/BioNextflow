/*
* braken subworkflow
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

process get_read_length {

  input:
  tuple val(pair_id), path(reads)

  output:
  stdout emit: out

  script:
  """
        if [ `echo ${reads} | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
        \$cat ${reads} | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", sum/100; exit} }'
  """
}

process bracken_build {

  maxForks 1

  tag { read_size }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  val(read_size)
  path(database)
  path(brackendb)

  output:
  path("database*"), emit: brackendb
  path("out${read_size}"), emit: bracken_out

  script:
  """
    if [ -f "${brackendb}/database.kraken" ]; then
      ln -s "${brackendb}/database.kraken" .
    else
      kraken2 --db=${database} --threads=${task.cpus} <( find -L ${database}/library \\( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \\) -exec cat {} + ) > database.kraken
    fi
    if [ -f "${brackendb}/database${read_size}mers.kmer_distrib}" ]; then
      ln -s "${brackendb}/database${read_size}mers.kraken" .
      ln -s "${brackendb}/database${read_size}mers.kmer_distrib" .
      touch out${read_size}
    else
      /usr/local/bracken/src/kmer2read_distr --seqid2taxid ${database}/seqid2taxid.map --taxonomy ${database}/taxonomy --kraken database.kraken --output database${read_size}mers.kraken -l ${read_size} -t ${task.cpus}
      python /usr/local/bracken/src/generate_kmer_distribution.py -i database${read_size}mers.kraken -o database${read_size}mers.kmer_distrib > out${read_size}
    fi
  """

}

process bracken {

  tag { pair_id }
  label (params.LABEL)
  container params.CONTAINER
  if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

  input:
  tuple val(pair_id), path(reads)
  path(brackendb)
  path("kraken2_${pair_id}.report")
  path("kraken2_${pair_id}.out")
  path(bracken_out)

  output:
  path("bracken_${pair_id}.*.report"), emit: report
  path("bracken_${pair_id}.*.out"), emit: output
  path("kraken2_${pair_id}_bracken_species.report"), emit: default_report

  script:
  """
  FIRST=\$(echo ${reads} | head -n1 | awk '{print \$1;}')
  if [ `echo \$FIRST | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
  READSIZE=\$(\$cat \$FIRST | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", sum/100; exit} } ')
  bracken -d ${brackendb} -i kraken2_${pair_id}.report -o bracken_${pair_id}.\${READSIZE}.report -r \$READSIZE -l S -t ${task.cpus} > bracken_${pair_id}.\${READSIZE}.out
  """

  }

  workflow BUILD {
      take:
      fastq
      database
      brackendb

      main:

      read_size = get_read_length(fastq).out.unique().map { it.trim().toInteger() }
      out = bracken_build(read_size, database, brackendb)

      emit:
      out.brackendb
      out.bracken_out

  }


  workflow RUN {
      take:
      fastq
      brackendb
      kraken2_report
      kraken2_outfile
      bracken_out

      main:
      out = bracken(fastq, brackendb.collect(), kraken2_report, kraken2_outfile, bracken_out.collect() )

      emit:
      out.report
      out.output
      out.default_report

  }

  workflow GET_VERSION {
      main:
  		getVersion()
      emit:
      getVersion.out
  }
