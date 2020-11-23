/*
*  lumpy module + samtools 
*/

params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = "bwa_out"
params.CONTAINER = "quay.io/biocontainers/lumpy-sv:0.3.0--py27hfbaaabd_6"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    lumpy 2>&1| grep Program | awk '{print "lumpy\t"\$4}'| sed s/\\)//g
    """
}

process getDiscordant {
    publishDir(params.OUTPUT, mode:'copy')

    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    publishDir("${params.OUTPUT}/lumpy_bam", mode:'copy')

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}.discordants.bam") 
 
    """    
    samtools view -@ ${task.cpus} -b -F 1294 ${pair_id}.bam > ${pair_id}.discordants.unsorted.bam
    samtools sort -@ ${task.cpus} -o ${pair_id}.discordants.bam ${pair_id}.discordants.unsorted.bam 
	rm ${pair_id}.discordants.unsorted.bam
    """    
}

process getSplitReads {
    publishDir(params.OUTPUT, mode:'copy')

    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    publishDir("${params.OUTPUT}/lumpy_bam", mode:'copy')

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}.splitters.bam") 
 
    """    
	samtools view -@ {task.cpus} -h ${reads} \
    	| extractSplitReads_BwaMem -i stdin \
    	| samtools view -@ {task.cpus} -Sb - > ${pair_id}.splitters.unsorted.bam
    samtools sort -@ ${task.cpus} -o ${pair_id}.splitters.bam ${pair_id}.splitters.unsorted.bam 
	rm ${pair_id}.splitters.unsorted.bam
    """    
}

process lumpy_express_single {
    publishDir("${params.OUTPUT}/lumpy_vcf", mode:'copy')
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER

    input:
    tuple val(pair_id), path(bam), path(discordant), path(split)

    output:
    tuple val(pair_id), path("${pair_id}.vcf") 
 
    """    
	lumpyexpress \
    	-B ${bam} \
    	-S ${split} \
    	-D ${discordant} \
    	-o ${pair_id}.vcf
    """    
}

workflow LUMPY_ALL_SINGLE {
    take: 
    input
    
    main:
		discordant = getDiscordant(input)
		split = getSplitReads(input)
		out = lumpy_express_single(input.join(discordant).join(split))
    emit:
    	out
}





workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
