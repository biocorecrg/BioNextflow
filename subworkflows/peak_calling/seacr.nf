/*
*  SEACR 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "biocorecrg/seacr:1.3"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo seacr 1.3
    """
}

process peakCall {
    label (params.LABEL)
    tag { comp_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(comp_id), path(sample), path(input)

    output:
    tuple val(comp_id), path("${comp_id}*.bed"), emit: peaks
    
	script:
    """
    SEACR.sh ${sample} ${input} ${params.EXTRAPARS} ${comp_id}
    """
}

process makeBedGraph {

    label (params.LABEL)
    tag { id }
    container params.CONTAINER

    input:
    tuple val(id), path(sample)
    file(genome)

    output:
    tuple val(id), path("${id}.fragments.bedgraph"), emit: bedgraph
    
	script:
    """
	bedtools bamtobed -bedpe -i ${sample} > ${id}.bed
	awk '\$1 != "." && \$1==\$4 && \$6-\$2 < 1000 {print \$0}' ${id}.bed > ${id}.clean.bed
	cut -f 1,2,6 ${id}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${id}.fragments.bed
	bedtools genomecov -bg -i ${id}.fragments.bed -g ${genome} > ${id}.fragments.bedgraph
	"""

}

workflow MAKE_GRAPHS {
    take: 
    sample
    genome
    
    main:
		makeBedGraph(sample, genome)
    emit:
    	bedgraph = makeBedGraph.out.bedgraph
}


workflow CALL {
    take: 
    comparisons
    
    main:
		peakCall(comparisons)
    emit:
    	peaks = peakCall.out.peaks
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
