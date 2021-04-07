/*
*  salmon module
* INPUT IS a channel with [ val(ID), [READL1, READL2, ...], optional [READR1, READR2, ...] ]
*/

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/salmon:1.2.1--hf69c8f4_0"
params.EXTRAPARS = ""


process index {
    label (params.LABEL)
    tag { reference }
    container params.CONTAINER

    input:
    path(reference)
    val(indexname)
    val(extrapars)

    output:
    path(indexname)
    
    """
    salmon index ${extrapars} -p ${task.cpus} -t ${reference} -i ${indexname}
    """
}

process mapSE {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER

    input:
    tuple val(pair_id), path(reads)
    path(index)
    val(extrapars)

    output:
    tuple val(pair_id), path("${pair_id}") 
    
	script:
    """
    salmon quant ${extrapars} --validateMappings --seqBias -l A --gcBias -p ${task.cpus} -i ${index} -r ${reads} -o ${pair_id}
    """
}

process mapPE {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER

    input:
    tuple val(pair_id), path(readsA), path(readsB)
    path(index)
    val(extrapars)

    output:
    tuple val(pair_id), path("${pair_id}") 
    
	script:
    """
    salmon quant ${extrapars} --validateMappings --seqBias -l A --gcBias -p ${task.cpus} -i ${index} -1 ${readsA} -2 ${readsB} -o ${pair_id}
    """
}

workflow MAPPER_SALMON {
    take: 
    action
    input
    index
    extrapars
    
    main:
		if(action == "index") {
			out = index(input, index, extrapars)
		} else if (action == "mapSE") {
			out =  mapSE(input, index, extrapars)
		} else if (action == "mapPE") {
			out =  mapPE(input, index, extrapars)
		} 
    emit:
    	out
}

