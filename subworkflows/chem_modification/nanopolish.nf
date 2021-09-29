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
    nanopolish eventalign -t ${task.cpus} --reads ${fastq} --bam ${bam} --genome ${reference} --samples --print-read-names --scale-events --samples 2>/dev/null | pigz -p ${task.cpus} > ${idsample}_${fast5_file}_event_align.tsv.gz
    """
}

/*
* CONCAT_EVENTS
*/

process concatenate_events {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:params.OUTPUTMODE , pattern: '*.gz') }

    container params.CONTAINER
    label (params.LABEL)
    tag "${idsample}" 
	
    input:
    tuple val(idsample), path("event_align_*") 
    
    output:
    tuple val(idsample), path("${idsample}_event_collapsed"), emit: cat_events
    tuple val(idsample), path("${idsample}_processed_perpos_median.tsv.gz"), emit: cat_medians


    script:
    """
	zcat event_align* | awk '!(/^contig/ && NR>1)' | tee   >(pigz -p ${task.cpus} -9 - > ${idsample}_combined.eventalign.tsv.gz) | NanopolishComp Eventalign_collapse -t ${task.cpus} -o ${idsample}_event_collapsed
	mean_per_pos.py -i ${idsample}_combined.eventalign.tsv.gz -o ${idsample} -s 500000
	pigz -p ${task.cpus} -9 ${idsample}_processed_perpos_median.tsv
    """
}


/*
* Estimate polyA tail size with nanopolish
*/
process polyAtail {
    container params.CONTAINER
    label (params.LABEL)
    tag "${sampleID}--${fast5}" 
  
	input:
	tuple val(sampleID), path(fast5), path(alignment), path(alnindex), path(fastq)
	path(reference)

	output:
	tuple val(sampleID), path("*.polya.estimation.tsv")

	script:
	def fast5_index=fast5.getSimpleName()
	"""
	#index reads
	nanopolish index -d ./ ${fastq}
	# polya length estimation
	nanopolish polya -r ${fastq} ${params.EXTRAPARS} -g ${reference} -t ${task.cpus} -b ${alignment} > ${sampleID}-${fast5_index}.polya.estimation.tsv
	"""
} 

process collect_polyA_results {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT,pattern: "*.polya.estimation.tsv", mode:params.OUTPUTMODE ) }
	tag { sampleID }  
	
	input:
	tuple val(sampleID), path("nanopol_*")
	
	output:
    tuple val(sampleID), path("${sampleID}.nanopol.len"), emit: filtered_est
    path("*.polya.estimation.tsv"), emit: polya_est

	script:
	"""
	cat nanopol_* | awk '!(NR>1 && /leader_start/)' | grep -v "READ_FAILED_LOAD"  >> ${sampleID}.polya.estimation.tsv	
	awk -F"\t" '{if (\$10=="PASS") print \$1"\t"\$9}' ${sampleID}.polya.estimation.tsv > ${sampleID}.nanopol.len
	"""

}




/*
*/

workflow POLYA_LEN {

    take: 
    fast5_files
    bams
    bais
    fastqs
    reference
   
    main:
    	datafiles = bams.join(bais).join(fastqs)
	   	polyAs = polyAtail(fast5_files.combine(datafiles, by: 0), reference)
		collect_polyA_results(polyAs.groupTuple())
	emit:
		filtered_est = collect_polyA_results.out.filtered_est
		polya_est = collect_polyA_results.out.polya_est
 }

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
    	datafiles = bams.join(bais).join(fastqs).join(summaries).join(indexes)
	   	aligned_events = eventalign(fast5_files.combine(datafiles, by: 0), reference)
		concat_events = concatenate_events(aligned_events.groupTuple()).cat_events
	emit:
		concat_events
 }


/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
