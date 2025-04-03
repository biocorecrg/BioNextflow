/*
* Sniffle subworkflows 
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
params.CONTAINER = "quay.io/biocontainers/sniffles:2.5.2--pyhdfd78af_0"

include { PossiblyUnzipGenome } from '../misc/misc.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    sniffles --version
    """
}

process sniffles {
    label (params.LABEL)
    tag "${idfile}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(idfile), path(bamfile), path(index)
	path(reference)  

    output:
    tuple val(idfile), path("${idfile}.vcf"), emit: vcf
    tuple val(idfile), path("${idfile}.snf"), emit: snf
    
    script:
    """
    sniffles ${params.EXTRAPARS} \
    --input ${bamfile} --vcf ${idfile}.vcf --snf ${idfile}.snf \
    --reference ${reference} --threads ${task.cpus} 
    """
}

process sniffles_multisample {
    label (params.LABEL)
    tag "${idfile}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(idfile), path(snvfiles)

    output:
    tuple val(idfile), path("${idfile}.vcf"), emit: vcf
    
    script:
	def snv_list = snvfiles.join(" ")
    """
    sniffles ${params.EXTRAPARS} \
    --input ${snv_list} --vcf ${idfile}.vcf --threads ${task.cpus} 
    """
}

workflow RUN {
    take: 
    input
    index
    reference
    
    main:
		sniffles(input.join(index), reference)
	emit:
    	vcf = sniffles.out.vcf
    	snf = sniffles.out.snf
	
}

workflow RUN_MULTI {
    take: 
    input
    
    main:
		sniffles_multisample(input)
	emit:
    	out = sniffles_multisample.out.vcf
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
