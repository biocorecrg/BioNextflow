/* 
 * Repository of functions about aligners
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class aligner {

	/* 
	 * Function for mapping SE or PE reads with STAR mapper. It reads both gzipped and plain fastq
 	 */
	
    static def mappingWithSTAR( seq_id, indexGenome, reads, cpus, extrapars="", debug="no") { 
        """
        	if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi
			if [ `echo ${reads} | grep ".gz"` ]; then gzip="--readFilesCommand zcat"
			else gzip=""
			fi
            	\$print STAR --genomeDir ${STARgenome} \
                     --readFilesIn ${reads} \
                     \$gzip \
                     --outSAMunmapped None \
                     --outSAMtype BAM SortedByCoordinate \
                     --runThreadN ${cpus} \
                     --quantMode GeneCounts \
                     ${extrapars} \
                     --outFileNamePrefix ${seq_id}

                \$print mkdir STAR_${seq_id}
                \$print mv ${seq_id}Aligned* STAR_${seq_id}/.
                \$print mv ${seq_id}SJ* STAR_${seq_id}/.
                \$print mv ${seq_id}ReadsPerGene* STAR_${seq_id}/.
                \$print mv ${seq_id}Log* STAR_${seq_id}/.   
        """
    }

}