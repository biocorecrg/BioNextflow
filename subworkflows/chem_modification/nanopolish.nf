/*
* NanoPolish
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = 	"biocorecrg/mopnanopolish:0.2"


/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		nanopolish --version | grep version      
    """
}


/*
*/

process index {

    container params.CONTAINER
    label (params.LABEL)
    tag "${idsample}" 
 		
    input:
    tuple val(idsample), path(fast5_folder), path(fastq), path(seqsummary)
    
    output:
    tuple val(idsample), path ("${fastq}.*") 

    script:
    
    """ 
      nanopolish index -s ${seqsummary} -d ${fast5_folder} ${fastq}
    """

}

/*
*/
process eventalign {
    container params.CONTAINER
    label (params.LABEL)
    tag "${idsample}--${fast5_file}" 
 	
    input:
    tuple val(idsample), path(fast5_file), path(bam), path(bai), path(fastq), path(seqsummary), path("*")
    path(reference)
    
    output:
    tuple val(idsample), path("*_event_align.tsv.gz")

    script:
    """ 
    mkdir ${idsample}/
    cd ${idsample}/; ln -s ../*.fast5 .; cd ../
    nanopolish eventalign -t ${task.cpus} --reads ${fastq} --bam ${bam} --genome ${reference} --samples --print-read-names --scale-events --samples | pigz -p ${task.cpus} > ${idsample}_${fast5_file}_event_align.tsv.gz
    """
}

/*
* FROM HERE


process concat_events {
//	publishDir outputNanopolish, pattern: "*_combined.eventalign.tsv.gz",  mode: 'copy'
//    publishDir outputNanopolish, pattern: "*_processed_perpos_median.tsv.gz", mode: 'copy'
    container params.CONTAINER
    label (params.LABEL)
    tag "${idsample}" 
	
    input:
    set idsample, file("event_align_*") 
    //from np_eventalign_folders.groupTuple() 
    
    output:
    tuple val(idsample), path("${idsample}_event_collapsed_align.tsv") 
    //into np_event_collapsed 
    path("${idsample}_combined.eventalign.tsv.gz") 
    path("${idsample}_processed_perpos_median.tsv.gz")


    script:
    """
	zcat event_align* | awk '!(/^contig/ && NR>1)' | tee   >(pigz -p ${task.cpus} -9 - > ${idsample}_combined.eventalign.tsv.gz) | NanopolishComp Eventalign_collapse -t ${task.cpus} -o ${idsample}_event_collapsed_align.tsv

	mean_per_pos.py -i ${idsample}_combined.eventalign.tsv.gz -o ${idsample} -s 500000
	pigz -p ${task.cpus} -9 ${idsample}_processed_perpos_median.tsv
         
    """
}

*/


/*
*/

workflow EVENTALIGN {

    take: 
    fast5_folders
    bams
    bais
    fastqs
    summaries
    reference
    
    main:
        indexes = index(fast5_folders.join(fastqs).join(summaries))
    	fast5_folders.map {
    		[it[0], file("${it[1]}/*.fast5")]
    	}.transpose().set{fast5_files}
    	datafiles = bams.join(bais).join(summaries).join(fastqs).join(indexes)
	   	//fast5_files.combine(datafiles, by: 0).view()
	   	out = eventalign(fast5_files.combine(datafiles, by: 0), reference)

	emit:
		out
 }


/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
