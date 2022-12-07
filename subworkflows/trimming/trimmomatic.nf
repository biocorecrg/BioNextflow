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
    tuple val (pair_id), path ("${pair_id}.trim_out.log"), emit: trim_log

    script:

    """
    trimmomatic PE \
    -threads ${task.cpus} -trimlog ${pair_id}.log \
    ${fastq} \
    ${pair_id}_1P.fq.gz \
    ${pair_id}_2P.fq.gz \
    ${pair_id}_1UP.fq.gz \
    ${pair_id}_2UP.fq.gz \
    ${params.EXTRAPARS} 2> ${pair_id}.trim_out.log
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
    tuple val (pair_id), path ("${pair_id}.trim_out.log"), emit: trim_log

    script:

    """
    trimmomatic SE \
    -threads ${task.cpus} -trimlog ${pair_id}.log \
    ${fastq} \
    ${pair_id}_P.fq.gz \
    ${params.EXTRAPARS} 2> ${pair_id}.trim_out.log
    """
}

process filterPEAdapter {
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }
   
    tag {pair_id }
    container params.CONTAINER
    label (params.LABEL)

    input:
    tuple val (pair_id), path(fastq), path(adapter)

    output:
    tuple val (pair_id), path ("${pair_id}_1P.fq.gz"), path ("${pair_id}_2P.fq.gz"), emit: trimmed_reads
    tuple val (pair_id), path ("${pair_id}.trim_out.log"), emit: trim_log

    script:

    """
    trimmomatic PE \
    -threads ${task.cpus} -trimlog ${pair_id}.log \
    ${fastq} \
    ${pair_id}_1P.fq.gz \
    ${pair_id}_2P.fq.gz \
    ${pair_id}_1UP.fq.gz \
    ${pair_id}_2UP.fq.gz \
    ${params.EXTRAPARS} 2> ${pair_id}.trim_out.log
    """
}

process filterSEAdapter {
    if (params.OUTPUT != "") {publishDir(params.OUTPUT, mode: 'copy') }
    
    tag {pair_id }
    container params.CONTAINER
    label (params.LABEL)

    input:
    tuple val (pair_id), path(fastq), path(adapter)

    output:
    tuple val (pair_id), path ("${pair_id}_P.fq.gz"), emit: trimmed_reads
    tuple val (pair_id), path ("${pair_id}.trim_out.log"), emit: trim_log

    script:

    """
    trimmomatic SE \
    -threads ${task.cpus} -trimlog ${pair_id}.log \
    ${fastq} \
    ${pair_id}_P.fq.gz \
    ${params.EXTRAPARS} 2> ${pair_id}.trim_out.log
    """
}



workflow FILTER {
    take: 
    fastq
    
    main:
    def sep_fastq = separateSEandPE(fastq)
        
    outpe = filterPE(sep_fastq.pe)
    outse = filterSE(sep_fastq.se)
    
    trimmed_reads = outpe.trimmed_reads.map{
    	[it[0], [it[1], it[2]] ] 
    }.mix(outse.trimmed_reads.map{
     	[it[0], [it[1]] ] 
    })

    emit:
        trimmed_reads = trimmed_reads 
        trim_log = outpe.trim_log.mix(outse.trim_log)       

}

workflow FILTERADAPTER {
    take: 
    fastq
    adapter
    
    main:
    def sep_fastq = separateSEandPE(fastq)
        
    outpe = filterPEAdapter(sep_fastq.pe.combine(adapter))
    outse = filterSEAdapter(sep_fastq.se.combine(adapter))
    trimmed_reads = outpe.trimmed_reads.map{
    	[it[0], [it[1], it[2]] ] 
    }.mix(outse.trimmed_reads.map{
     	[it[0], [it[1]] ] 
    })

    emit:
        trimmed_reads = trimmed_reads
        trim_log = outpe.trim_log.mix(outse.trim_log)       

}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 


