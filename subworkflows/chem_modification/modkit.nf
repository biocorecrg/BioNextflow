/*
* Epinano
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = "quay.io/biocontainers/ont-modkit:0.4.1--h5c23e0d_0"

include { FAIDX as SAMTOOLS_FAIDX } from "../misc/samtools" addParams(EXTRAPARS: "", OUTPUT: "")


/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		later
    """
}


/*
*/

process pileup {

    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }
    label (params.LABEL)
    tag "${id}" 
 	
    input:
    tuple val(id), path(bam), path(index)
    
    output:
    path("${id}_pileup.bed.gz")
    
    script:

	"""
        modkit pileup -t ${task.cpus}  ${bam} ${id}_pileup.bed
        gzip ${id}_pileup.bed
	"""
}



/*
*/

workflow PILEUP {

    take: 
    bam
    bai

    main:
    	out = pileup(bam.join(bai))
 		 	        
}

/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
