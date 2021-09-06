/*
* trimmomatic module
*/

include { separateSEandPE } from '../global_functions.nf'

params.CONTAINER = "quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4"
params.OUTPUT = "trimmomatic_output"

params.LABEL = ""
params.EXTRAPARS = ""


process getVersion {
    container params.CONTAINER

    output:
    stdout emit: out

    shell:
    """
    echo trimmomatic v0.32
    """
}

process filterPE {
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }
   
    tag {pair_id }
    container params.CONTAINER
    label (params.LABEL)

    input:
    tuple val (pair_id), path(fastq)

    output:
    tuple val (pair_id), path ("${pair_id}_1P.fq.gz"), path ("${pair_id}_2P.fq.gz"), emit: trimmed_reads
    tuple val (pair_id), path ("${pair_id}.log"), emit: trim_log

    script:

    """
    trimmomatic PE \
    -threads ${task.cpus} -trimlog ${pair_id}.log \
    ${fastq} \
    ${pair_id}_1P.fq.gz \
    ${pair_id}_2P.fq.gz \
    ${pair_id}_1UP.fq.gz \
    ${pair_id}_2UP.fq.gz \
    ${params.EXTRAPARS}
    """
}

process filterSE {
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }
    
    tag {pair_id }
    container params.CONTAINER
    label (params.LABEL)

    input:
    tuple val (pair_id), path(fastq)

    output:
    tuple val (pair_id), path ("${pair_id}_P.fq.gz"), emit: trimmed_reads
    tuple val (pair_id), path ("${pair_id}.log"), emit: trim_log

    script:

    """
    trimmomatic SE \
    -threads ${task.cpus} -trimlog ${pair_id}.log \
    ${fastq} \
    ${pair_id}_P.fq.gz \
    ${params.EXTRAPARS}
    """
}


workflow FILTER {
    take: 
    fastq
    
    main:
    def sep_fastq = separateSEandPE(fastq)
        
    outpe = filterPE(sep_fastq.pe)
    outse = filterSE(sep_fastq.se)

    emit:
        trimmed_reads = outpe.trimmed_reads.mix(outse.trimmed_reads)
        trim_log = outpe.trim_log.mix(outse.trim_log)       

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 


