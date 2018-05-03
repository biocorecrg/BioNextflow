/* 
 * Repository of functions about aligners
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class NGSaligner {

	/*
	 * Indexing a genome with Bowtie2 mapper. It reads both gzipped and plain fasta
 	 */
	
    static def indexWithBowtie2( genome_file, indexname="bowtie2genome", cpus, extrapars="", debug="no") { 
 
        """			    
    	if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi	

     	if [ `echo ${genome_file} | grep ".gz"` ]; then 
			\$print zcat ${genome_file} > `basename ${genome_file} .gz`
        	\$print bowtie2-build --threads ${cpus} `basename ${genome_file} .gz` ${indexname}
        	\$print rm `basename ${genome_file} .gz`
		else \$print bowtie2-build --threads ${cpus} ${genome_file} ${indexname}
		fi
        """
    }

	/*
	 * Mapping SE and PE reads with Bowtie2. Reads can be both gzipped and plain fastq
	*/ 	     
    static def mappingSEWithBowtie2(reads, indexGenome, alnfile, cpus, extrapars="", debug="no") { 
        """ 
    	if [ `echo ${debug} == "debug"` ]; then 
    	echo "bowtie2 --non-deterministic -x ${indexGenome} -U ${reads} -p ${cpus}" '| samtools view -Sb -@ ' ${cpus} '- >' ${alnfile}; 
    	else 
    	bowtie2 --non-deterministic -x ${indexGenome} -U ${reads} -p ${cpus} | samtools view -Sb -@ ${cpus} - > ${alnfile}
    	fi
    		
        """
	}
 	
	/* 
	 * Mapping SE or PE reads with STAR mapper. It reads both gzipped and plain fastq

	
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
 	 */
}