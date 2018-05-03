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



}
