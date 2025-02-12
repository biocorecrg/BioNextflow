/*
* MultiQC subworkflows 
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/quay.io/biocontainers/multiqc:1.27.1--pyhdfd78af_0"


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
	path("*_report.html"), emit: report
	path("*_data"), emit: data
	
	
    script:
    """
		multiqc ${params.EXTRAPARS} .
    """
}

process makeReportWithConfig {
    label (params.LABEL)
    
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    path(input)
    path(config)
	
    output:
	path("*_report.html"), emit: report
	path("*_data"), emit: data
	
	
    script:
    """
		multiqc ${params.EXTRAPARS} -c ${config} .
    """
}

process makeReportID {
    label (params.LABEL)
    tag { id }
    
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(input)
	
    output:
	tuple val(id), path(id)
	
    script:
    """
		multiqc ${params.EXTRAPARS} -o ${id} --filename ${id}_multiqc .
    """
}

workflow REPORT {
    take: 
    input
    
    main:
		makeReport(input)
	emit:
		out	= makeReport.out.report
		data = makeReport.out.data
}

workflow REPORT_WITH_CONFIG {
    take: 
    input
    config
    
    main:
		makeReportWithConfig(input, config)
	emit:
		out	= makeReportWithConfig.out.report
		data = makeReportWithConfig.out.data
}




workflow REPORT_ID {
    take: 
    input
    
    main:
		out = makeReportID(input)
	emit:
		out	
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
