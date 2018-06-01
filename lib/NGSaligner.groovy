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
	 * Indexing a transcriptome with Salmon mapper. It reads both gzipped and plain fasta
 	 */
	
    static def indexWithSalmon( transcript_file, indexname="transcript.index", kmer, cpus, extrapars="", debug="no") { 
 
        """	
		if [ `echo ${transcript_file} | grep ".gz"` ]; then 
			zcat ${transcript_file} > `basename ${transcript_file} .gz`;
   		    salmon index -t `basename ${transcript_file} .gz` -i ${indexname} --type quasi -k ${kmer} ${extrapars} -p ${cpus};
   		    rm `basename ${transcript_file} .gz`;
		else salmon index -t ${transcript_file} -i ${indexname} --type quasi -k ${kmer} ${extrapars} -p ${cpus}
		fi
        """
    }

    /*
	 * Mapping to a transcriptome index with Salmon mapper. 
 	 */
	
    static def mapPEWithSalmon( transcript_index, readsA, readsB, output, libtype="ISF", cpus, extrapars="", debug="no") { 
 
        """	               
		if [ `echo ${readsA} | grep ".gz"` ]; then 
			 salmon quant -i ${transcript_index} --gcBias -l ${libtype} -1 <(gunzip -c ${readsA}) -2 <(gunzip -c ${readsB}) ${extrapars} -o ${output}
		else  
			salmon quant -i ${transcript_index} --gcBias -l ${libtype} -1 ${readsA} -2 ${readsB} ${extrapars} -o ${output}
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
	 * Indexing a genome with STAR mapper. It reads both gzipped and plain fasta
	*/ 
	
    static def indexWithSTAR( genome_file, outfolder, outprefix, annotation, readsize, cpus, extrapars="") { 
        """
		mkdir ${outfolder}
     	if [ `echo ${genome_file} | grep ".gz"` ]; then 
			zcat ${genome_file} > `basename ${genome_file} .gz`
			STAR --runMode genomeGenerate --genomeDir ${outfolder} --runThreadN ${cpus} \
			--genomeFastaFiles `basename ${genome_file} .gz` --sjdbGTFfile ${annotation} \
			--sjdbOverhang ${readsize} --outFileNamePrefix ${outprefix};
			rm `basename ${genome_file} .gz`
		else 
			STAR --runMode genomeGenerate --genomeDir ${outfolder} --runThreadN ${cpus} \
			--genomeFastaFiles ${genome_file} --sjdbGTFfile ${annotation} \
			--sjdbOverhang ${readsize} --outFileNamePrefix ${outprefix}
		fi
		"""
    }
    
    /*
	 * Mapping SE and PE reads with START. Reads can be both gzipped and plain fastq
	*/ 	     
    static def mappingWithSTAR(pair_id, genome, reads, cpus, extrapars="") { 
    """
	if [ `echo "${reads}"| cut -f 1 -d " " | grep ".gz"` ]; then gzipped=" --readFilesCommand zcat "; else gzipped=""; fi
		STAR --genomeDir ${genome} \
				 --readFilesIn ${reads} \
				 \$gzipped \
				 --outSAMunmapped None \
				 --outSAMtype BAM SortedByCoordinate \
				 --runThreadN ${cpus} \
				 --quantMode GeneCounts \
				 --outFileNamePrefix ${pair_id}
			 
			mkdir STAR_${pair_id}
			mv ${pair_id}Aligned* STAR_${pair_id}/.
			mv ${pair_id}SJ* STAR_${pair_id}/.
			mv ${pair_id}ReadsPerGene* STAR_${pair_id}/.
			mv ${pair_id}Log* STAR_${pair_id}/.   
    """
	}
    
    
    /*
	 * Indexing bam alignments with samtools
 	 */
	
    static def indexAlnWithSamtools( bamfile) { 
 
        """
        samtools index ${bamfile}
        """
    }

}