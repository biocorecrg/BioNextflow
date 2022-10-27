/*
*  feature count module
*  This workflow allows to make count on input data
*  It needs input aln_data (id, bam and index)
*/

params.LABEL = ""
params.CONTAINER = "quay.io/biocontainers/subread:2.0.1--h7132678_2"
params.OUTPUT = ""
params.EXTRAPARS = ""

include { unzipNamedPipe } from '../global_functions.nf'



process featureCounts {
    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern:'*.counts', mode:'copy') }
   
    input:
    tuple val(id), path(bamfile), path(intervals)

    output:
    tuple val(id), path("${id}.counts"), emit: counts
    tuple val(id), path("${id}.log"), emit: logs
  
	script:
	if ("${intervals}" =~ /.bed$/ ) {
		precmd = "awk \'OFS=\"\t\" {print \$1\"-\"\$2\"-\"\$3, \$1, \$2, \$3, \".\"} \'  ${intervals} >  ${intervals}.saf; "
		interv = "${intervals}.saf"
		extrapars = " -F SAF "

	} else {
	    precmd = ""
		interv = "${intervals}"
		extrapars = ""
	}
	"""
		${precmd} \
	    featureCounts ${extrapars} ${params.EXTRAPARS} -a ${interv} -T ${task.cpus} -o ${id}.counts ${bamfile} 2> ${id}.log
	"""

}

workflow COUNT {
    take: 
	data
    
    main:    
    out = featureCounts(data)
    
    emit:
	counts = out.counts
	logs = out.logs

}

workflow MULTI_COUNT {
    take: 
    intervals
    bams
    
    main: 
    data = intervals.combine(bams.map{[it]}).map{ [ it[0], it[2], it[1]]}    
    out = featureCounts(data)
    
    emit:
	counts = out.counts
	logs = out.logs

}


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo "featureCounts -version"
    """
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
