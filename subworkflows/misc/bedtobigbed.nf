/*
*  samtools 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/ucsc-bedtobigbed@sha256:f27a6df27a91eb6ae64cbc18359c242b6642b0832e1b30737686b66b44157a66"
params.OUTPUTMODE = "copy"


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    bedToBigBed 2>&1| grep Conv | awk -F "-" '{print \$1}'
    """
}


process bedToBigBed {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    file(chromsizes)
    tuple val(pair_id), path(bed)

    output:
    tuple val(pair_id), path("${pair_id}.bb") 
    
	script:
    """    
    bedToBigBed ${bed} ${chromsizes} ${pair_id}.bb
    """
}


workflow BED2BIGBED {
    take: 
    genome
    bed
    
    main:
		out = bedToBigBed(genome, bed)
    emit:
    	out
}




workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
