/*
*  bwa module + samblaster + samtools 
*/

include { BWA_INDEX as BWA_INDEX } from "./bwa"

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "bwa_sv_out"
params.CONTAINER = "quay.io/biocontainers/mulled-v2-8a9a988fff4785176b70ce7d14ff00adccf8a5b8:aeac8200e5c50c5acf4dd14792fd8453255af835-0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    bwa 2>&1| grep Version | awk '{print "bwa "\$0}';
    samblaster -h 2>&1 | grep Version;
    samtools --version | grep samtools
    """
}

process mapPE {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER

    input:
    tuple val(pair_id), path(reads)
    path(indexes)

    output:
    tuple val(pair_id), path("${pair_id}.bam") 
    
	script:
    def indexname = indexes[0].baseName

    """    
    bwa mem -R "@RG\\tID:id\\tSM:sample\\tLB:lib" -t ${task.cpus} ${indexname} ${reads} \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view  -@ ${task.cpus} -S -b - > ${pair_id}.bam;
    """
}

process mapSE {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER

    input:
    tuple val(pair_id), path(reads)
    path(indexes)

    output:
    tuple val(pair_id), path("${pair_id}.bam") 
    
	script:
    def indexname = indexes[0].baseName

    """    
    bwa mem -R "@RG\\tID:id\\tSM:sample\\tLB:lib" -t ${task.cpus} ${indexname} ${reads} \
    | samblaster --excludeDups --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view  -@ ${task.cpus} -S -b - > ${pair_id}.bam;
    """
}


workflow SAMBLASTER_MAP_PE {
    take: 
    input
    indexes
    
    main:
		out = mapPE(input, indexes)
    emit:
    	out
}

workflow SAMBLASTER_MAP_SE {
    take: 
    input
    indexes
    
    main:
		out = mapSE(input, indexes)
    emit:
    	out
}


workflow SAMBLASTER_ALL_PE {
    take: 
    reference
    input
    
    main:
		index = BWA_INDEX(reference)
		out = SAMBLASTER_MAP_PE(input, index)
    emit:
    	out
}

workflow SAMBLASTER_ALL_SE {
    take: 
    reference
    input
    
    main:
		index = BWA_INDEX(reference)
		out = SAMBLASTER_MAP_SE(input, index)
    emit:
    	out
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
