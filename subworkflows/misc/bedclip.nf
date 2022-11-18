/*
*  samtools 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/ucsc-bedclip:377--h0b8a92a_2"

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


process bedClip {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    file(chromsizes)
    tuple val(pair_id), path(bed)

    output:
    tuple val(pair_id), path("${pair_id}_clip.bed")
   
    script:
    """
    bedClip ${bed} ${chromsizes} ${pair_id}_clip.bed
    """
}


workflow BEDCLIP {
    take: 
    genome
    bed
    
    main:
		out = bedClip(genome, bed)
    emit:
    	out
}




workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
