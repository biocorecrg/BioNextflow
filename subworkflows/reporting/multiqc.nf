/*
* MultiQC subworkflows 
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = "multiqc_out"
params.CONTAINER = "quay.io/biocontainers/multiqc:1.10.1--pyhdfd78af_1"


include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    multiqc --version
    """
}

process makeReport {
    label (params.LABEL)
    
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    path(input)
	
    output:
	path("multiqc_report.html")
	
    script:
    """
		multiqc .
    """
}


workflow REPORT {
    take: 
    input
    
    main:
		out = makeReport(input)
	emit:
		out	
}



workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
