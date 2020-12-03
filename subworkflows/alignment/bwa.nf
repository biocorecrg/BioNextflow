/*
*  bwa module + samtools 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "bwa_out"
params.CONTAINER = "quay.io/biocontainers/mulled-v2-8a9a988fff4785176b70ce7d14ff00adccf8a5b8:aeac8200e5c50c5acf4dd14792fd8453255af835-0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    bwa 2>&1| grep Version | awk '{print "bwa "\$0}';
    samtools --version | grep samtools
    """
}


process index {
    label (params.LABEL)
    tag { reference }
    container params.CONTAINER

    input:
    path(reference)

    output:
    path("${reference}.*")
    
    """
    bwa index ${params.EXTRAPARS} ${reference}
    """
}

process map {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    publishDir(params.OUTPUT, mode:'copy')

    input:
    tuple val(pair_id), path(reads)
    path(indexes)

    output:
    tuple val(pair_id), path("${pair_id}.bam") 
    
	script:
    def indexname = indexes[0].baseName

    """    
    bwa mem -t ${task.cpus} ${indexname} ${reads} | samtools view -@ ${task.cpus} -Sb > ${pair_id}.bam
    """
}

workflow BWA_MAP {
    take: 
    input
    indexes
    
    main:
		out = map(input, indexes)
    emit:
    	out
}

workflow BWA_INDEX {
    take: 
    reference
    
    main:
		ref_file = file(reference)
		if( !ref_file.exists() ) exit 1, "Missing ${reference} file!"
		out = index(reference)

    emit:
    	out
}

workflow BWA_ALL {
    take: 
    reference
    input
    
    main:
		index = BWA_INDEX(reference)
		out = BWA_MAP(input, index)
    emit:
    	out
}



workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
