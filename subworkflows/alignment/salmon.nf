/*
*  salmon module
* INPUT IS a channel with [ val(ID), [READL1, READL2, ...], optional [READR1, READR2, ...] ]
*/


include { separateSEandPE } from '../global_functions.nf'

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/salmon:1.2.1--hf69c8f4_0"
params.EXTRAPARS = ""

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
    tag { reference }
    container params.CONTAINER

    input:
    path(reference)
    val(indexname)

    output:
    path(indexname)
    
    """
    salmon index ${params.EXTRAPARS} -p ${task.cpus} -t ${reference} -i ${indexname}
    """
}

process mapSE {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy'}

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
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy'}

    input:
    tuple val(pair_id), path(reads)
    path(index)

    output:
    tuple val(pair_id), path("${pair_id}") 
    
    script:
    """
    salmon quant ${params.EXTRAPARS} --validateMappings --seqBias -l A --gcBias -p ${task.cpus} -i ${index} -1 ${reads[0]} -2 ${reads[1]} -o ${pair_id}
    """
}


workflow INDEX {
    take: 
    reference
    
    main:
	ref_file = file(reference)
	if( !ref_file.exists() ) exit 1, "Missing ${reference} file!"
	def refname = ref_file.simpleName
        out = index(reference, refname)
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
    fastq
    
    main:        
    index = INDEX(reference)
    out = MAP(index, fastq)


    emit:
        out     

}

