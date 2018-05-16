/* 
 * Repository for biological functions
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class Misc {

 	 
	/* 
	 * Function for running fastQC on input samples
 	 */	
    static def samtoolSort(bamfile, sortedfile, cpus="1", debug="no") {

    """
    	if [ `echo ${debug} == "debug"` ]; then 
    	echo samtools sort -@ ${cpus} ${bamfile} '>' ${sortedfile}; else 	
		samtools sort -@ ${cpus} ${bamfile} > ${sortedfile};
		fi
    """
	}
	
	/* 
	 * Function for running fastQC on input samples
 	 */
    static def getTranscriptsFromGTF(genome_file, annotation, output="transcript.fa", debug="no") {

    """
    	if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi	
     	if [ `echo ${genome_file} | grep ".gz"` ]; then 
			\$print zcat ${genome_file} > `basename ${genome_file} .gz`
			\$print gffread -g `basename ${genome_file} .gz` -w ${output} ${annotation}
        	\$print rm `basename ${genome_file} .gz`
        else \$print gffread -g ${genome_file} -w ${output} ${annotation}
		fi
    """
	}
	




}
