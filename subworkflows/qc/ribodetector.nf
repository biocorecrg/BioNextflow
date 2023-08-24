/*
*  QC module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.CONTAINER = "quay.io/biocontainers/ribodetector:0.2.7--pyhdfd78af_0"
params.OUTPUT = ""

include { CALC_AVG_READSIZE } from "../misc/misc.nf"
include { DOWNSAMPLE_PAIRS } from "../misc/misc.nf"
include { unzipCmd }  from "../misc/misc.nf"


process riboDetector {
    tag "${id}"
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(id), path(fastqs)
    val(readsize)

    output:
	tuple val(id), path ("${id}_rna_num.txt")

    script:

    def unzipA      = unzipCmd(fastqs[0])
    def file_nameA  = unzipA[0]
    def cmd_unzipA  = unzipA[1]
    def cmd_cleanA  = unzipA[2]
    def unzipB      = unzipCmd(fastqs[1])
    def file_nameB  = unzipB[0]
    def cmd_unzipB  = unzipB[1]
    def cmd_cleanB  = unzipB[2]
    """
	${cmd_unzipA}
	${cmd_unzipB}
     ribodetector_cpu  -t ${task.cpus} -l ${readsize} -o no1.fq.gz no2.fq.gz -e rrna -r rna1.fq rna2.fq -i ${file_nameA} ${file_nameB} 
	 awk -v id=${id} '{num++}END{print id" "num/4}' rna1.fq > ${id}_rna_num.txt
      ${cmd_cleanA} 
      ${cmd_cleanB}
    """
}

process makeMultiQCReport {

    container params.CONTAINER
   

    input:
 	path(reports)
    val(downsize)
    
    output:
	path ("ribo_stats_mqc.txt")

    script:

	"""
echo '# id: ribodetector
# plot_type: bargraph
# section_name: Ribosome contamination
# description: % of ribosomal reads on ${downsize}) 
Filename	reads' > ribo_stats_mqc.txt;
	awk -v tot=${downsize} '{print \$1"\t"\$2/tot*100}' ${reports} >> ribo_stats_mqc.txt
	"""
	
	
}

workflow QC {
    take: 
    fastqs
    downsize
    
    main:
    avg_size = CALC_AVG_READSIZE(fastqs, "illumina")
    down_reads = DOWNSAMPLE_PAIRS(fastqs, downsize)
    out = riboDetector(down_reads, avg_size).map{it[1]}.collect()
	mqc = makeMultiQCReport(out, downsize)
	 
    emit:
    out = out
    mqc = mqc

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    ribodetector_cpu -v
    """
}



workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
