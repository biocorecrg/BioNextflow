/*
*  Homer 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "biocorecrg/chipanno:0.3"

include { unzipCmd } from '../global_functions.nf'

process addGeneNames {
    label (params.LABEL)
    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(anno)
    path(annofile)

    output:
    tuple val(id), path("*_g.anno"), emit: g_anno
    tuple val(id), path("*.stats"), emit:  stats
    
	script:
    """
    addGeneNameToHomer.pl -g ${annofile} -a ${anno} -o ${id}_g.anno
    awk -F "\\t" '{split(\$8,a," "); if (a[1]!="Annotation") {print a[1]}}' ${id}_g.anno |sort |uniq -c| awk -v file=${id} 'BEGIN{print "Features\\t"file}{print \$2"\\t"\$1}' > ${id}.anno.stats 
    """
}

process makeAnnoReport {
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    path(anno_vals)

    output:
 	file("peak_anno_stats_mqc.txt")
    
	script:
    """
    makePeakAnnoReport.pl -ext .anno.stats -o peak_anno_stats_mqc.txt
    """
}


workflow MAKE_REPORT {
    take: 
    anno_peaks
    
    main:
		out = makeAnnoReport(anno_peaks)
    emit:
    	out 
}


workflow ADD_GENE_NAMES {
    take: 
    anno_peaks
    annofile
    
    main:
		addGeneNames(anno_peaks, annofile)
    emit:
    	out = addGeneNames.out.g_anno
    	stats = addGeneNames.out.stats
}


