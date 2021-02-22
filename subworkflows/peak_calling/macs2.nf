/*
*  Macs2 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "moaims_out"
params.CONTAINER = "quay.io/biocontainers/macs2:2.2.7.1--py37h516909a_0"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    macs2 --version
    """
}


process peakCall {
    label (params.LABEL)
    tag { comp_id }
    container params.CONTAINER
    publishDir(params.OUTPUT, mode:'copy')

    input:
    tuple val(comp_id), path(sample), path(input)
    val(gsize)

    output:
    tuple val(comp_id), path("${comp_id}_peaks.narrowPeak"), emit: narrowPeaks
    tuple val(comp_id), path("${comp_id}_peaks.xls"), emit: xlsPeaks
    tuple val(comp_id), path("${comp_id}_summits.bed"), emit: bedSummits
    
	script:

    """
    macs2 callpeak ${params.EXTRAPARS} -t ${sample} -c ${input} -g ${gsize} -n ${comp_id} --fix-bimodal --call-summits
    """
}


workflow MACS2_CALL {
    take: 
    comparisons
    gsize
    
    main:
		peakCall(comparisons, gsize)
    emit:
    	narrowPeaks = peakCall.out.narrowPeaks
    	xlsPeaks = peakCall.out.xlsPeaks
    	bedSummits = peakCall.out.bedSummits
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
