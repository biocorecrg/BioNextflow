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
	
	/* 
	 * Function for making coverage profiles
	 * WARNING there is no distinction between strands and it will be replaced in near future with something better

    static def makeAlnProfiles(bamfile, readsize, genome_file, output="profile.bw") {

	"""
	ratio=`samtools idxstats $1| grep -v '*' | awk -v readsize=$2 '{sum+=\$3}END{print 1000000000/(sum*readsize)}'`;
	echo $ratio > ratio.txt
	bedtools genomecov -bg -split -ibam $1 -g $3 -scale $ratio > `basename $1`.bg
	bedSort `basename $1`.bg `basename $1`.bg
	bedGraphToBigWig `basename $1`.bg $3 $4
	rm `basename $1`.bg
	"""
	}
	*/
}
