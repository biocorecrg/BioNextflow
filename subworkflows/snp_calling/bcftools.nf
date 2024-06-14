/*
* ClairS subworkflows 
*
* The parameters are: 
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for star
*	OUTPUT for storing the final sub-workflow output 
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/bcftools:1.20--h8b25389_0"


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    bcftools --version| head -n 1
    """
}

process view {
    label (params.LABEL)
    tag "${id}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.filt.vcf.gz")
    
    script:
    """
		bcftools view --threads ${task.cpus} -o ${id}.filt.vcf.gz  -O z ${params.EXTRAPARS} ${vcf}  
    """
}



workflow VIEW {

    take: 
    vcf
    
    main:
        out = view(vcf)
        
	emit:
    	out
	
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
