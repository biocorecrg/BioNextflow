/* 
 * Class for misc functions (to be splitted?) 
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class Misc {

	/*
	 * Properties definition
	 */
	
     String input = ''
     String output = ''
     String mode = ''
     Integer number = 0
     Integer read_size = 0
     String memory = 0
     String chr_size_file = ''
     Integer cpus = 1
     String extrapars = ''

	/* 
	 *  Sorting bam files with samtools
 	 */	
    def public st_sortBam() {

    """
		samtools sort -@ ${this.cpus} ${this.input} > ${this.output};
    """
	}

    /*
	 * Indexing bam alignments with samtools
 	 */
	
 	def public st_indexBam() { 
 
        """
        samtools index ${this.input}
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
	 * Function for making coverage profiles. It needs samtools, bedtools, bedSort and bedGraphToBigWig
	 * WARNING there is no distinction between strands and it will be replaced in near future with something better
     */

    def public makeAlnProfiles() {

	"""
	ratio=`samtools idxstats ${this.input}| grep -v '*' | awk -v readsize=${this.read_size} '{sum+=\$3}END{print 1000000000/(sum*readsize)}'`;
	echo \$ratio > ratio.txt
	bedtools genomecov -bg -split -ibam ${this.input} -g ${this.chr_size_file} -scale \$ratio > `basename ${this.input}`.bg
	bedSort `basename ${this.input}`.bg `basename ${this.input}`.bg
	bedGraphToBigWig `basename ${this.input}`.bg ${this.chr_size_file} ${this.output}
	rm `basename ${this.input}`.bg
	"""
	}

}
