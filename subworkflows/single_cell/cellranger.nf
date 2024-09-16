/*
*  cellranger module
*/

params.LABEL = ""
params.CONTAINER = "biocorecrg/debian-perlbrew:buster"
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
    cellranger -V
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

    cellranger mkgtf ${params.EXTRAPARSINDEX} ${anno_name} anno.filtered.gtf 
    cellranger mkref \
 		 --genome=${indexname} \
 		 --fasta=${genome_name} \
  		 --genes=anno.filtered.gtf \
  		 --nthreads=${task.cpus} \
  		 --memgb=${task.memory.toGiga()}

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
    ln -s ${pairs[0]} ./${pair_id}_S1_L001_R1_001.fastq.gz
    ln -s ${pairs[1]} ./${pair_id}_S1_L001_R2_001.fastq.gz
    
	cellranger count ${params.EXTRAPARS} --id=${pair_id} \
                   --transcriptome=${index} \
                   --fastqs=./ \
                   --sample=${pair_id} \
                   --localcores=${task.cpus} \
                   --localmem=${task.memory.toGiga()}   
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


