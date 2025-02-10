params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/slow5tools:1.3.0--h56e2c18_0"


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    samtools --version | grep samtools
    """
}

process pod5_2_blow5 {

    label (params.LABEL)
    tag "${ idfile }"
    container params.CONTAINER
		
	input:
	tuple val(idfile), path(pod5)

	output:
	tuple val(idfile), path("${idfile}.blow5")

	
	script:
	"""
		blue-crab p2s -t ${task.cpus} ${pod5} -o ${idfile}.blow5
	"""
}


process blow5_merge {

    label (params.LABEL)
    tag "${ idfile }"
    container params.CONTAINER
		
	input:
	tuple val(idfile), path(pod5s)

	output:
	tuple val(idfile), path("${idfile}.blow5")

	
	script:
	"""
		slow5tools merge -t ${task.cpus} ./ -o ${idfile}.blow5
	"""
}




workflow CONVERT {

    take: 
    input_pod5    
    
    main:	   
   	  out = pod5_2_blow5(input_pod5)            

	emit:
		out
 
 }
 
 workflow MERGE {

    take: 
    input_pod5    
    
    main:	   
   	  out = blow5_merge(input_pod5)            

	emit:
		out
 
 }
 
 


