/* 
 * Repository for biological functions
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class BiologicalFunctions {

	/* 
	 * Function for mapping PE reads with STAR mapper. It reads both gzipped and plain fastq
 	 */
	
    static String mappingPairsWithSTAR( pair_id, STARgenome, reads, cpus, debug="no") { 
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
                     --outFileNamePrefix ${pair_id}

                \$print mkdir STAR_${pair_id}
                \$print mv ${pair_id}Aligned* STAR_${pair_id}/.
                \$print mv ${pair_id}SJ* STAR_${pair_id}/.
                \$print mv ${pair_id}ReadsPerGene* STAR_${pair_id}/.
                \$print mv ${pair_id}Log* STAR_${pair_id}/.   
        """
        .stripIndent()
    }

	/* 
	 * Function for trimming read pairs with Skewer
 	 */
	
    static def trimPairsWithSkewer( pair_id, reads, minread, cpus, debug="no") { 
		"""
			skewer -t ${cpus} -l ${minread} -n -u -o ${pair_id} -z ${reads}
    	"""
	}




}