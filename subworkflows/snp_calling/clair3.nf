/*
* Clair3 subworkflows 
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
params.CONTAINER = "hkubal/clair3:v1.0.10"

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

process clair3 {
    label (params.LABEL)
    tag "${comp_id}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(comp_id), path(bam) , path(bai), path(reference), path(refai), path("*")

    output:
    tuple val(comp_id), path(comp_id), emit: folder
    tuple val(comp_id), path("${comp_id}/*.vcf.gz"), emit: vcf
    tuple val(comp_id), path("${comp_id}/*.bam"), path("${comp_id}/*.bai"), emit: phased_bam
    
    script:
    """

    /opt/bin/run_clair3.sh \
	--bam_fn ${bam} \
    --ref_fn ${reference} \
    --output ${comp_id} \
    --threads ${task.cpus} ${params.EXTRAPARS} 
    rm -fr ${comp_id}/tmp
    """

}



workflow RUN {

    take: 
    bams
    reference
    reffai
    model
    
    main:
        ref_genome = CHECK_FASTA(reference)
        models = model.map {
            [ it ] 
        }
        mydata = bams.combine(ref_genome).combine(reffai).combine(models)
       
        out = clair3(mydata)
        
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
