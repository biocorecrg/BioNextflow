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
params.CONTAINER = "quay.io/biocontainers/snpeff:5.2--hdfd78af_1"


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    snpEff -version
    """
}

process snpeff_ann {
    label (params.LABEL)
    tag "${id}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.vcf.gz') }
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.genes.txt') }

    
    input:
    tuple val(id), path(vcf), val(genomeid), path(spneffdata), path(config)

    output:
    tuple val(id), path("${id}.ann.vcf.gz"), emit: annotation
    tuple val(id), path("snpEff_summary.html"), emit: html_summary
    tuple val(id), path("${id}.snpEff_summary.csv"), emit: csv_summary
    tuple val(id), path("${id}.snpEff_summary.genes.txt"), emit: genes
    
    script:
    """
    snpEff -Xmx${task.memory.giga}g ann -csvStats ${id}.snpEff_summary.csv -dataDir ./${spneffdata} ${genomeid} ${vcf} | gzip -c > ${id}.ann.vcf.gz
    """
}


workflow ANN {

    take: 
    vcf
    genomeid
    genomedata
    config
    
    main:
        genomeid_ch = channel.value(genomeid)
        out = snpeff_ann(vcf.combine(genomeid_ch).combine(genomedata).combine(config))
        
	emit:
    	annotation = out.annotation
    	html_summary = out.html_summary
    	csv_summary = out.csv_summary
    	genes = out.genes
	
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
