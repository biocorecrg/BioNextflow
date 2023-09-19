/*
* Trust4 subworkflows 
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/trust4:1.0.12--h43eeafb_0"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    run-trust4 2>&1| head -n 1| cut -d " " -f 1,2
    """
}

process run_se {
    label (params.LABEL)
    tag "${id}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(reads), path(rec_genes), path(reference)

    output:
    tuple val(id), path("${id}*"), emit: results 
         
	script:
    """    
     run-trust4 -f ${rec_genes} -t ${task.cpus} --ref ${reference} -u ${reads} -o ${id} ${params.EXTRAPARS}
    """
}

workflow RUN {
    take: 
    reference
    rec_genes
    reads
    
    main:
		out = run_se(reads.combine(rec_genes).combine(reference))
	emit:
    	out
 	
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
