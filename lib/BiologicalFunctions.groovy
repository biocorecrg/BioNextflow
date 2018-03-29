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
	
    static String mappingPairsWithSTAR( pair_id, STARgenome, reads, cpus, extrapars="", debug="no") { 
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
	
    static def trimPairsWithSkewer( pair_id, reads, minread, cpus, extrapars="", debug="no") { 
	"""
		if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi
		\$print	skewer ${extrapars} -t ${cpus} -l ${minread} -n -u -o ${pair_id} -z ${reads}
    	"""
        .stripIndent()
     }


	/* 
	 * Function for mapping single ends with Bowtie2. It returns a bam file ${read_id}.bam (conversion via samtools)
 	 */
	
    static def mappingSingleEndsWithBowtie2(read_id, single_reads, index, cpus, extrapars="", debug="no") { 
	"""
		if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi
		\$print	"bowtie2 --non-deterministic ${extrapars} -x ${index} -U ${single_reads} -p ${cpus} | samtools view -Sb -@ ${cpus} - > ${read_id}.bam"
    """
        .stripIndent()
     }

	/* 
	 * Function for indexing a genome via Bowtie1. It returns an index composed by a number of files genome_bowtie1*  
	 * It reads both gzipped and plain fastq
 	 */
    static def indexingGenomeWithBowtie1(genomeFile, single_reads, index, cpus, extrapars="", debug="no") { 
	"""
		if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi
		if [ `echo ${genomeFile} | grep 'gz'` ]; then \$print "zcat ${genomeFile} > genome_bowtie1.fa"; else \$print ln -s ${genomeFile} genome_bowtie1.fa; fi
		\$print bowtie-build ${extrapars} --threads ${cpus} genome_bowtie1.fa genome_bowtie1
		\$print rm genome_bowtie1.fa
    """
        .stripIndent()
     }	




}
