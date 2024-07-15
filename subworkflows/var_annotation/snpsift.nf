/*
* snpeff subworkflows 
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
params.CONTAINER = "quay.io/biocontainers/snpsift:5.2--hdfd78af_0"


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    snpsift -version
    """
}

process snpsift_ann {
    label (params.LABEL)
    tag "${id}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.vcf.gz') }

    
    input:
    tuple val(id), path(vcf), path(vcf_ref), path(vcf_ref_idx)

    output:
    tuple val(id), path("*.vcf.gz")
    
    script:
    def name_ref = vcf_ref.simpleName
    """
    SnpSift annotate ${vcf_ref} ${vcf} > ${id}.on.${name_ref}.vcf
    gzip ${id}.on.${name_ref}.vcf
    """
}


workflow ANN {

    take: 
    vcf
    vcf_ref
    vcf_ref_idx
    
    main:
        out = snpsift_ann(vcf.combine(vcf_ref).combine(vcf_ref_idx))
        
	emit:
    	out
	
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
