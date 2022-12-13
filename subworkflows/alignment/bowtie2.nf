/*
* bowtie subworkflows (samtools embedded for converting the output) 
* The parameters are: 
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output 
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.LABELINDEX = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    bowtie2 --version | head -n 1 | awk '{print "bowtie2 "\$3}' 
    samtools --version | grep samtools
    """
}


process index {
    label (params.LABEL)
    label (params.LABELINDEX)
    
    tag { "${reference}" }
    container params.CONTAINER

    input:
    path(reference)

    output:
    path("${reference}.*.bt2")
    
    """
    bowtie2-build --threads ${task.cpus} -f ${reference} ${reference}
    """
}

process map {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads)
    path(indexes)

    output:
    tuple val(pair_id), path("${pair_id}.bam") 
    
	script:
    def indexname = "${indexes[0]}".replaceAll(".1.bt2", "")
    def cmd = "-U ${reads[0]}"
    def cmd2 = ""
	if (reads[1]) {
		cmd = "-1 ${reads[0]} -2 ${reads[1]}"
	}
	if (reads[2] && reads[3]) {
		cmd2 = "-U ${reads[2]},${reads[3]}"
	}

    """    
    bowtie2 -x ${indexname} -p ${task.cpus} ${cmd} ${cmd2} ${params.EXTRAPARS} | samtools view -@ ${task.cpus} -Sb > ${pair_id}.bam
    """
}

workflow MAP {
    take: 
    input
    indexes
    
    main:
		out = map(input, indexes)
    emit:
    	out
}

workflow INDEX {
    take: 
    reference
    
    main:
		ref_file = file(reference)
		if( !ref_file.exists() ) exit 1, "Missing ${reference} file!"
		out = index(reference)

    emit:
    	out
}

workflow ALL {
    take: 
    reference
    input
    
    main:
		index = INDEX(reference)
		out = MAP(input, index)
    emit:
    	out
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
