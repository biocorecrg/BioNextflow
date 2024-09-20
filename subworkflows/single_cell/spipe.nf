/*
*  spipe module
*/

params.LABEL = ""
params.CONTAINER = "biocorecrg/spipe:1.3.1"
params.EXTRAPARS = ""
params.EXTRAPARSINDEX = ""
params.OUTPUT = ""

include { separateSEandPE } from '../global_functions.nf'
include { unzipCmd } from '../global_functions.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out

    shell:
    """
    split-pipe --version
    """
}

process index {
    label (params.LABEL)
    tag { "${indexname}" }
    container params.CONTAINER

    input:
    path(annotation)
    path(genome)
    val(indexname)

    output:
    path(indexname)

    script:

    def unzip_genome = unzipCmd(genome)
    def cmd_genome = unzip_genome[1]
    def genome_name = unzip_genome[0]
    def clean_genome = unzip_genome[2]

    def unzip_annot = unzipCmd(annotation)
    def cmd_annot = unzip_annot[1]
    def anno_name = unzip_annot[0]
    def clean_anno = unzip_annot[2]
 

    """
    ${cmd_genome}
    ${cmd_annot}

    split-pipe \
	   --mode mkref {params.EXTRAPARSINDEX} \
	   --genome_name ${indexname}  \
	   --fasta ${genome_name} \
       --genes ${anno_name} \
       --nthreads ${task.cpus} \
	   --output_dir ./indexname

    ${clean_genome}
    ${clean_anno}    
    
    """
}


process count {
    label (params.LABEL)
    tag { "${pair_id}" }

    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy')}

    input:
    tuple val(pair_id), path(pairs), path(index)

    output:
    tuple val(pair_id), path("${pair_id}")

    script:    
    """
    split-pipe \
       --mode all ${params.EXTRAPARS} \
       --genome_dir ${index} \
       --fq1 ${pairs[0]} \
       --fq2 ${pairs[1]} \
       --output_dir ${pair_id} \
    
    """
}



workflow INDEX {
    take:
    annotation
    genome

    main:
		ref_file = file(genome)
		if( !ref_file.exists() ) exit 1, "Missing ${genome} file!"
		def refname = ref_file.simpleName
        out = index(annotation, genome, refname)

    emit:
    	out
}

workflow COUNT {
    take:
    index
    fastq

    main:
    out = count(fastq.combine(index))


    emit:
        out 

}


workflow ALL {
    take:
    annotation
    genome
    fastq

    main:
    index = INDEX(annotation, genome)
    out = COUNT(index, fastq)


    emit:
    	index
    	out

}


