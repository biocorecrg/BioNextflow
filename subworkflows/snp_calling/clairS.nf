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
params.CONTAINER = "hkubal/clairs:v0.2.0"

include { PossiblyUnzipGenome } from '../misc/misc.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    clairs --version
    """
}

process clairS {
    label (params.LABEL)
    tag "${idfile}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple path(bamcancer), path(bamctrl), path(reference)

    output:
    //tuple val(idfile), path("${idfile}.snf"), emit: snf, optional true
    
    script:
    """

    /opt/bin/run_clairs \
        --tumor_bam_fn ${bamcancer} \
        --normal_bam_fn ${bamctrl} \
        --ref_fn ${reference} \
        --threads ${task.cpus} ${params.EXTRAPARS} \
        --output_dir ./ \
        --conda_prefix /opt/conda/envs/clairs

    """
}


workflow RUN {

    take: 
    bams
    reference
    
    main:
        ref_genome = PossiblyUnzipGenome(reference)
        clairS(bams.combine(ref_genome))
        
	//emit:
    	//vcf = clairS.out.vcf
    	//snf = clairS.out.snf
	
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
