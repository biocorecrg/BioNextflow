/*
* Misc subworkflows 
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "biocontainers/biocontainers:v1.2.0_cv1"

include { unzipCmd } from '../global_functions.nf'


process makePeakReport {
    label (params.LABEL)
    
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    path(input)
	
    output:
	path("peak_stats_mqc.txt")
	
    script:
    """
echo '# id: peakcall
# plot_type: bargraph
# section_name: Peak calling statistics
# description: Number of peaks
Filename	peaks' > peak_stats_mqc.txt;
	for i in `ls ${input}`; do wc -l \$i | awk '{print \$2"\\t"\$1}' | sed s/_peaks\\.bed//g >> peak_stats_mqc.txt; done
    """
}


workflow PEAK_REPORT {
    take: 
    input
    
    main:
		out = makePeakReport(input)
	emit:
		out	
}




