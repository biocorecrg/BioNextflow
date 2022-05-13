/*
*  deeptools
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/deeptools:3.5.1--py_0"
params.OUTPUTMODE = "copy"

include { unzipCmd } from '../global_functions.nf'


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    deeptools --version
    """
}


process BamCoverageChipSeq {
    label (params.LABEL)

    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(id), path(bam), path(bai)
    val(effgsize)

    output:
    tuple val(id), path("${id}.bw") 
    
	script:
    """    
	bamCoverage --bam ${bam} -o ${id}.bw \
   		--binSize 10 \
    	--normalizeUsing RPGC \
    	--effectiveGenomeSize ${effgsize} \
		${params.EXTRAPARS} \
		-p ${task.cpus}
    """
}



workflow BAMCOV_CHIP {
    take: 
    bam
    gsize
    
    main:
		out = BamCoverageChipSeq(bam, gsize.first())
    emit:
    	out
}




workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
