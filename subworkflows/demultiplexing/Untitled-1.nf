 /* process extracting_demultiplexed_fastq {
		label 'basecall_cpus'
   	    tag {"${demultiplexer}"}  
				
		input:
    	set idfile, file(demux), file(fastq) from demux_for_fastq_extraction.join(fastq_files_for_demultiplexing)
        
		output:
		set idfile, file ("*.fastq.gz") into fastq_for_filtering

		script:
		"""
            extract_sequence_from_fastq.py ${demux} ${fastq}
			for i in *.fastq; do gzip \$i; done
 		"""
	} 
	
process extracting_demultiplexed_fast5 {

		label 'basecall_cpus'
   	    tag { demultiplexer }  
		publishDir outputFast5,  mode: 'copy'
				
		input:
    	file("demux_*") from demux_for_fast5_extraction.collect()
    	file("*") from fast5_files_for_demultiplexing.collect()

		output:
        file("*")
        
		script:
		"""
		cat demux_* | grep -v ReadID >> dem.files
		awk '{print \$2 > \$3".list" }' dem.files
		for i in *.list; do mkdir `basename \$i .list`; fast5_subset --input ./ --save_path `basename \$i .list`/ --read_id_list \$i --batch_size 4000 -t ${task.cpus}; done 
 		rm *.list
 		rm * / filename_mapping.txt
 		rm dem.files 
 		"""
	} 

 */

 process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		#guppy_basecaller --version       
    """
}