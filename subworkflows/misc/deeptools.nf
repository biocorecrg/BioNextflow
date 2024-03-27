/*
*  deeptools
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/deeptools:3.5.1--py_0"
params.OUTPUTMODE = "copy"

include { unzipCmd } from '../global_functions.nf'


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    deeptools --version
    """
}


process BamCoverageChipSeq {
    label (params.LABEL)

    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(id), path(bam), path(bai)
    val(effgsize)

    output:
    tuple val(id), path("${id}.bw") 
    
	script:
    """    
	bamCoverage --bam ${bam} -o ${id}.bw \
   		--binSize 10 \
    	--normalizeUsing RPGC \
    	--effectiveGenomeSize ${effgsize} \
		${params.EXTRAPARS} \
		-p ${task.cpus}
    """
}

process BamCoverageChipSeqScale {
    label (params.LABEL)

    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(id), path(bam), path(bai), val(scale)
    val(effgsize)

    output:
    tuple val(id), path("${id}.bw") 
    
	script:
    """    
	bamCoverage --bam ${bam} -o ${id}.bw \
   		--binSize 10 \
    	--normalizeUsing RPGC \
    	--scaleFactor ${scale} \
    	--effectiveGenomeSize ${effgsize} \
		${params.EXTRAPARS} \
		-p ${task.cpus}
    """
}

process calcFingerPrints {
    label (params.LABEL)

    //tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE, pattern: "*.pdf") }

    input:
    tuple path(bam), path(bai)

    output:
    path("fingerprints.pdf"), emit: pdf 
    path("fingerprint.qcmetrics.txt"), emit: metrics 
    path("fingerprint.count.txt"), emit: counts 
    
	script:
    """    
	plotFingerprint \
	-b ${bam} \
	--smartLabels \
	--minMappingQuality 30 --skipZeros \
	--numberOfSamples 50000 \
	-T "Fingerprints of samples"  \
    --plotFile fingerprints.pdf \
    --outQualityMetrics fingerprint.qcmetrics.txt \
    --outRawCounts fingerprint.count.txt \
	${params.EXTRAPARS} \
	-p ${task.cpus}
    """
}


process computeMatrixForGenes {
    label (params.LABEL)

    container params.CONTAINER

    input:
    path(bigwigs)
    path(annotation_gtf)

    output:
    path("matrix.gz") 
    
	script:
    """    
	computeMatrix scale-regions -S ${bigwigs} \
	-R ${annotation_gtf}  \
	--smartLabels  \
	${params.EXTRAPARS}   \
	--beforeRegionStartLength 3000    \
	--regionBodyLength 5000  \
	--afterRegionStartLength 3000  \
	--numberOfProcessors ${task.cpus}   \
	--skipZeros -o matrix.gz
    """
}

process computeSingleMatrixForGenes {
    label (params.LABEL)
    tag { id }

    container params.CONTAINER

    input:
    tuple val(id), path(bigwigs)
    path(annotation_gtf)

    output:
    tuple val(id), path("${id}_genes_matrix.gz") 
    
	script:
    """    
	computeMatrix scale-regions -S ${bigwigs} \
	-R ${annotation_gtf}  \
	--smartLabels  \
	${params.EXTRAPARS}   \
	--beforeRegionStartLength 3000    \
	--regionBodyLength 5000  \
	--afterRegionStartLength 3000  \
	--numberOfProcessors ${task.cpus}   \
	--skipZeros -o ${id}_genes_matrix.gz
    """
}


process computeMatrixForTSS {
    label (params.LABEL)

    container params.CONTAINER

    input:
    path(bigwigs)
    path(annotation_gtf)

    output:
    path("matrix.gz") 
    
	script:
    """    
	computeMatrix reference-point -S ${bigwigs} \
	-R ${annotation_gtf} --referencePoint TSS \
	--smartLabels  \
	${params.EXTRAPARS}   \
	--beforeRegionStartLength 3000    \
	--afterRegionStartLength 3000  \
	--numberOfProcessors ${task.cpus}   \
	--skipZeros -o matrix.gz
    """
}

process computeSingleMatrixForTSS {
    label (params.LABEL)
    tag { id }

    container params.CONTAINER

    input:
    tuple val(id), path(bigwigs)
    path(annotation_gtf)

    output:
    tuple val(id), path("${id}_TSS_matrix.gz") 
    
	script:
    """    
	computeMatrix reference-point -S ${bigwigs} \
	-R ${annotation_gtf} --referencePoint TSS \
	--smartLabels  \
	${params.EXTRAPARS}   \
	--beforeRegionStartLength 3000    \
	--afterRegionStartLength 3000  \
	--numberOfProcessors ${task.cpus}   \
	--skipZeros -o ${id}_TSS_matrix.gz
    """
}


process plotTSSprofile {
    label (params.LABEL)

    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    path(matrix)

    output:
    path("enrichment_TSS.pdf") 
    
	script:
    """    
	plotHeatmap -m ${matrix}  \
	-out enrichment_TSS.pdf  --heatmapWidth 20\
	--colorMap jet  --missingDataColor "#FFF6EB" --heatmapHeight 15
    """
}

process plotSingleTSSprofile {
    label (params.LABEL)
    tag { id }

    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(id), path(matrix)

    output:
    path("${id}_TSS.pdf") 
    
	script:
    """    
	plotHeatmap -m ${matrix}  \
	-out ${id}_TSS.pdf  --heatmapWidth 10\
	--colorMap jet  --missingDataColor "#FFF6EB" --heatmapHeight 10
    """
}

process plotGeneProfile {
    label (params.LABEL)

    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    path(matrix)

    output:
    path("enrichment_Genes.pdf") 
    
	script:
    """    
	plotHeatmap -m ${matrix}  \
	-out enrichment_Genes.pdf  --heatmapWidth 20\
    --perGroup -T "Read enrichment in gene body"    
    """
}

process plotSingleGeneProfile {
    label (params.LABEL)
    tag { id }

    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(id), path(matrix)

    output:
    tuple val(id), path("${id}_genes.pdf") 
    
	script:
    """    
	plotHeatmap -m ${matrix}  \
	-out ${id}_genes.pdf  --heatmapWidth 20\
    --perGroup -T "Read enrichment in gene body"    
    """
}

process BamCoverage {
    label (params.LABEL)

    tag { id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE) }

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id), path("${id}.bw") 
    
	script:
    """    
	bamCoverage --bam ${bam} -o ${id}.bw ${params.EXTRAPARS} -p ${task.cpus}
    """
}



workflow CALC_FINGERPRINGS {

    take: 
    bams
    bai
    
    main:
		out = calcFingerPrints(bams.mix(bai).collate(2))
    
    emit:
    	pdf = out.pdf 
		metrics = out.metrics 
    	counts = out.counts 
}

workflow BAMCOV_CHIP_SCALE {
    take: 
    bam
    gsize
    
    main:
		out = BamCoverageChipSeqScale(bam, gsize)
    emit:
    	out
}


workflow BAMCOV_CHIP {
    take: 
    bam
    gsize
    
    main:
		out = BamCoverageChipSeq(bam, gsize)
    emit:
    	out
}

workflow BAMCOVERAGE {
    take: 
    bam
    
    main:
		out = BamCoverage(bam)
    emit:
    	out
}

workflow PLOT_COV_TSS {
    take: 
    bigwigs
    annotation_gtf
    
    main:
		matx = computeMatrixForTSS(bigwigs, annotation_gtf)
		out = plotTSSprofile(matx)
    emit:
    	out
}

workflow PLOT_SCOV_TSS {
    take: 
    bigwigs
    annotation_gtf
    
    main:
		matx = computeSingleMatrixForTSS(bigwigs, annotation_gtf)
		out = plotSingleTSSprofile(matx)
    emit:
    	out
}

workflow PLOT_COV_GENES {
    take: 
    bigwigs
    annotation_gtf
    
    main:
		matx = computeMatrixForGenes(bigwigs, annotation_gtf)
		out = plotGeneProfile(matx)
    emit:
    	out
}

workflow PLOT_SCOV_GENES {
    take: 
    bigwigs
    annotation_gtf
    
    main:
		matx = computeSingleMatrixForGenes(bigwigs, annotation_gtf)
		out = plotSingleGeneProfile(matx)
    emit:
    	out
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
