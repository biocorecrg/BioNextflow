/*
*  cellranger_multi module
*/

params.LABEL = ""
params.CONTAINER = "quay.io/nf-core/cellranger:9.0.1"
params.EXTRAPARS = ""
params.EXTRAPARSINDEX = ""
params.OUTPUT = ""

include { separateSEandPE } from '../global_functions.nf'
include { unzipCmd } from '../global_functions.nf'
include { INDEX } from './cellranger.nf'

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out

    shell:
    """
    cellranger -V
    """
}



process count {
    label (params.LABEL)
    tag { "${pair_id}" }

    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy')}

    input:
    tuple val(pair_id), path(pairs), val(barcodes),  path(index)

    output:
    tuple val(pair_id), path("${pair_id}")

    script:
	def list_ids = barcodes.join("\n")  

    """
    ln -s ${pairs[0]} ./${pair_id}_S1_L001_R1_001.fastq.gz
    ln -s ${pairs[1]} ./${pair_id}_S1_L001_R2_001.fastq.gz
    echo  "[gene-expression]" > multi_config.csv
    echo  "reference,\$PWD/${index}" >> multi_config.csv
    echo  "create-bam,false" >> multi_config.csv
    echo  "" >> multi_config.csv
    echo  "[libraries]" >> multi_config.csv
    echo  "fastq_id,fastqs,feature_types" >> multi_config.csv
    echo  "${pair_id},\$PWD,Gene Expression" >> multi_config.csv
    echo  "" >> multi_config.csv
    echo  "[samples]" >> multi_config.csv
    echo  "sample_id,ocm_barcode_ids" >> multi_config.csv
    echo  "${list_ids}" >> multi_config.csv

    cellranger multi --id=${pair_id} --csv=multi_config.csv  \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()}     

    """
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


