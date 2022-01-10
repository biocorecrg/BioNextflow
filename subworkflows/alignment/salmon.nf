/*
*  salmon module
* INPUT IS a channel with [ val(ID), [READL1, READL2, ...], optional [READR1, READR2, ...] ]
*/

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/salmon:1.5.1--h84f40af_0"
params.EXTRAPARS = ""
params.EXTRAPARSINDEX = ""
params.OUTPUT = ""

include { separateSEandPE } from '../global_functions.nf'
include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out

    shell:
    """
    salmon --version
    """
}

process index {
    label (params.LABEL)
    tag { "${reference}" }
    container params.CONTAINER

    input:
    path(reference)
    path(genome)
    val(indexname)

    output:
    path("${indexname}")

    script:
    def unzip_ref = unzipCmd(reference)
    def cmd_ref = unzip_ref[1]
    def ref_name = unzip_ref[0]
    def clean_ref = unzip_ref[2]
    def unzip_genome = unzipCmd(genome)
    def cmd_genome = unzip_genome[1]
    def genome_name = unzip_genome[0]
    def clean_genome = unzip_genome[2]
    """
    ${cmd_ref}
    ${cmd_genome}

    grep "^>" ${genome_name} | cut -d " " -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat ${ref_name} ${genome_name} >> gentrome.fa
    gzip gentrome.fa
    salmon index ${params.EXTRAPARSINDEX} -p ${task.cpus} -d decoys.txt -t gentrome.fa.gz -i ${indexname}

    ${clean_ref}
    ${clean_genome}
    rm gentrome.fa.gz
    """
}

process mapSE {
    label (params.LABEL)
    tag { "${pair_id}" }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy')}

    input:
    tuple val(pair_id), path(reads)
    path(index)

    output:
    tuple val(pair_id), path("${pair_id}")

    script:
    """
    salmon quant ${params.EXTRAPARS} --validateMappings --seqBias -l A --gcBias -p ${task.cpus} -i ${index} -r ${reads} -o ${pair_id}
    """
}

process mapPE {
    label (params.LABEL)
    tag { "${pair_id}" }

    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy')}

    input:
    tuple val(pair_id), path(pairs)
    path(index)

    output:
    tuple val(pair_id), path("${pair_id}")

    script:
    def readsA = pairs[0]
    def readsB = pairs[1]
    """
    salmon quant ${params.EXTRAPARS} --validateMappings --seqBias -l A --gcBias -p ${task.cpus} -i ${index} -1 ${readsA} -2 ${readsB} -o ${pair_id}
    """
}


workflow INDEX {
    take:
    reference
    genome

    main:
	ref_file = file(reference)
	if( !ref_file.exists() ) exit 1, "Missing ${reference} file!"
	def refname = ref_file.simpleName
        out = index(reference, genome, refname)
    emit:
    	out
}

workflow MAP {
    take:
    index
    fastq

    main:
    def sep_fastq = separateSEandPE(fastq)
    outpe = mapSE(sep_fastq.se, index)
    outse = mapPE(sep_fastq.pe, index)


    emit:
        out = outpe.mix(outse)

}

workflow ALL {
    take:
    reference
    genome
    fastq

    main:
    index = INDEX(reference, genome)
    outm = MAP(index, fastq)


    emit:
    	index
    	outm

}
