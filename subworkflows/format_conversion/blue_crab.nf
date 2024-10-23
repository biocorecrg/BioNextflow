params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.LABEL = ""
params.CONTAINER = "biocorecrg/blue-crab:0.2.0"


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





workflow CONVERT {

    take: 
    input_pod5    
    
    main:	   
   	  out = pod5_2_blow5(input_pod5)            

	emit:
		out
 
 }
 
 


