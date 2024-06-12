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

include { CHECK_FASTA } from '../misc/misc.nf'

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
    tag "${comp_id}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(comp_id), path(bamcancer), path(bamctrl), path(baicancer), path(baictrl), path(reference), path(refai)

    output:
    tuple val(comp_id), path(comp_id), emit: folder
    tuple val(comp_id), path("${comp_id}/output.vcf.gz"), emit: vcf
    
    script:
    """

    /opt/bin/run_clairs \
        --tumor_bam_fn ${bamcancer} \
        --normal_bam_fn ${bamctrl} \
        --ref_fn ${reference} \
        --threads ${task.cpus} ${params.EXTRAPARS} \
        --output_dir ${comp_id} \
        --conda_prefix /opt/conda/envs/clairs

    """
// --remove_intermediate_dir
}



workflow RUN {

    take: 
    bams
    reference
    reffai
    
    main:
        ref_genome = CHECK_FASTA(reference)
        out = clairS(bams.combine(ref_genome).combine(reffai))
        
	emit:
    	vcf = out.vcf
    	folder = out.folder
	
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
