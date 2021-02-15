/*
* STAR subworkflows 
* The accessible subworkflows are:
* GET_VERSION that emits the version of bwa and samtools as stdout
* STAR_INDEX that takes:
*	a channel with an optionally gzipped fasta file
*   it emits a list of files as index
* STAR_MAP that takes:
*	a channel list with index files as produced by STAR_INDEX
*	a channel containing one or two (gzipped) fastq files 
*	it emits a channel containing a tuple of id, bam file
* STAR_ALL (STAR_INDEX + STAR_MAP) that takes:
*   a channel containing one or two (gzipped) fastq files
*	a channel with an optionally gzipped fasta file
*	a channel with an optionally gzipped gtf file (or "" if not available)
*	a channel with an overhang value (or "" if has to be calculated)
*   it emits a channel containing a tuple of id, bam file
*
* The parameters are: 
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output 
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = "bwa_out"
params.CONTAINER = "quay.io/biocontainers/star:2.7.7a--0"

include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    STAR --version
    """
}

process calcOverhang {
    label (params.LABEL)
    
    tag { pair_id }
    container params.CONTAINER

    input:
    tuple val(pair_id), path(reads)
	
    output:
    stdout emit: readsize

    script:
    def first_pair = reads[0]
    """
    if [ `echo ${first_pair} | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
        \$cat ${first_pair} | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", (sum/100)-1; exit} }'
    """
}


process indexWithAnno {
    label (params.LABEL)
    tag { "indexing ${reference} with ${annotation} and overhang=${overhang}"  }
    container params.CONTAINER

    input:
    val(overhang)
    val(indexname)
    path(reference)
	path(annotation)  

    output:
    path("${indexname}")
    
    script:
    def unzip_ref = unzipCmd(reference)
    def cmd_ref = unzip_ref[1]
    def ref_name = unzip_ref[0]
    def clean_ref = unzip_ref[2]
    def unzip_anno = unzipCmd(annotation)
    def cmd_anno = unzip_anno[1]
    def anno_name = unzip_anno[0]
    def clean_anno = unzip_anno[2]
    """
		${cmd_ref} 
		${cmd_anno}
        mkdir ${indexname}
        STAR --runMode genomeGenerate --genomeDir ${indexname} \
        	--sjdbOverhang ${overhang} --sjdbGTFfile ${anno_name} \
            --runThreadN ${task.cpus} \
            --genomeFastaFiles ${ref_name} \
            --outFileNamePrefix ${indexname}
        ${clean_ref}
        ${clean_anno}
    """
}

process indexNoAnno {
    label (params.LABEL)
    tag { reference }
    container params.CONTAINER

    input:
    val(indexname)
    path(reference)

    output:
    path("${indexname}")
    
    script:
    def unzip_ref = unzipCmd(reference)
    def cmd_ref = unzip_ref[1]
    def ref_name = unzip_ref[0]
    def clean_ref = unzip_ref[2]
    """
		${cmd_ref} 
        mkdir ${indexname}
        STAR --runMode genomeGenerate --genomeDir ${indexname} \
            --runThreadN ${task.cpus} \
            --genomeFastaFiles ${ref_name} \
            --outFileNamePrefix ${indexname}
        ${clean_ref}
    """
}


process map {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    publishDir(params.OUTPUT, mode:'copy', pattern: '*.bam')

    input:
    path(indexes)
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}*.bam"), emit: bams 
    tuple val(pair_id), path("${pair_id}PerGene*"), emit: quants optional true
    tuple val(pair_id), path("${pair_id}Log*"), emit: logs 
    tuple val(pair_id), path("${pair_id}SJ*"),  emit: junctions
         
	script:
    """    
        if [ `echo "${reads}"| cut -f 1 -d " " | grep ".gz"` ]; then gzipped=" --readFilesCommand zcat "; else gzipped=""; fi
            STAR --genomeDir ${indexes} \
                 --readFilesIn ${reads} \
                  \$gzipped \
                  --outSAMunmapped None \
                  --outSAMtype BAM SortedByCoordinate \
                  --runThreadN ${task.cpus} \
                  --outFileNamePrefix ${pair_id} \
                  ${params.EXTRAPARS}
   """
}

workflow STAR_MAP {
    take: 
    input
    indexes
    
    main:
		map(indexes, input)
	emit:
    	bams = map.out.bams
    	logs = map.out.logs
    	quants = map.out.quants
    	junctions = map.out.junctions
	
}

workflow STAR_INDEX {
    take: 
	reference
    annotation
	overhang    
    
    main:
		ref_file = file(reference)
		anno_file = file(annotation)
		if( !ref_file.exists() ) exit 1, "Missing ${reference} file!"
		if( !anno_file.exists() ) exit 1, "Missing ${anno_file} file!"
		def refname = ref_file.simpleName
		out = indexWithAnno(overhang, refname, reference, annotation)
    emit:
    	out
}

workflow STAR_INDEX_NOANNO {
    take: 
	reference
    
    main:
		ref_file = file(reference)
		if( !ref_file.exists() ) exit 1, "Missing ${reference} file!"
		def refname = ref_file.simpleName
		out = indexNoAnno(refname, reference)
    emit:
    	out
}

workflow STAR_ALL {

    take: 
    input
    reference
    annotation
    overhang
    
    main:
		if (overhang == "") {    
			overhang_val = calcOverhang(input.first())
		} else {
			overhang_val = overhang
		}
		if (annotation != "") {
			index = STAR_INDEX(reference, annotation, overhang_val)
		} else {
			index = STAR_INDEX_NOANNO(reference)
		}
		STAR_MAP(input, index)
	emit:
		bams = STAR_MAP.out.bams
		logs = STAR_MAP.out.logs
		quants = STAR_MAP.out.quants
		junctions = STAR_MAP.out.junctions
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
