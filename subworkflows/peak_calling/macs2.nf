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
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(comp_id), path(sample), path(input)
    val(gsize)

    output:
    tuple val(comp_id), path("${comp_id}_peaks.narrowPeak"), optional: true, emit: narrowPeaks
    tuple val(comp_id), path("${comp_id}_peaks.xls"), emit: xlsPeaks
    tuple val(comp_id), path("${comp_id}_summits.bed"), optional: true,  emit: bedSummits
    tuple val(comp_id), path("${comp_id}_peaks.gappedPeak"), optional: true, emit: gappedPeaks
    tuple val(comp_id), path("${comp_id}_peaks.broadPeak"), optional: true, emit: broadPeaks
    
	script:

    """
    macs2 callpeak ${params.EXTRAPARS} -t ${sample} -c ${input} -g ${gsize} -n ${comp_id} --fix-bimodal
    """
}


workflow CALL {
    take: 
    comparisons
    gsize
    
    main:
		peakCall(comparisons, gsize)
    emit:
    	narrowPeaks = peakCall.out.narrowPeaks
    	broadPeaks = peakCall.out.broadPeaks
    	xlsPeaks = peakCall.out.xlsPeaks
    	bedSummits = peakCall.out.bedSummits
    	gappedPeaks = peakCall.out.gappedPeaks
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
