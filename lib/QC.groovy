/* 
 * Repository for biological functions
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class QC {

	/* 
	 * Function for mapping PE reads with STAR mapper. It reads both gzipped and plain fastq
	*/
	
    static def qualimapRNAseq ( bamfile, annotation_file, outfolder="QUALIMAP", strand="strand-specific-reverse", memory="2G", single="YES", extrapars="", debug="no") { 

    """
       if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi
       if [ `echo ${single} == "YES"` ]; then pe_mode="-pe"; else pe_mode=""; fi
      
	   \$print unset DISPLAY
       \$print mkdir tmp
       \$print export JAVA_OPTS="-Djava.awt.headless=true -Xmx${memory} -Djava.io.tmpdir=\$PWD/tmp"    
       \$print qualimap rnaseq \${pe_mode} --java-mem-size=${memory} -bam ${bamfile} -gtf ${annotation_file} -outdir ${outfolder} -p ${strand} ${extrapars}
       \$print rm -fr tmp
    """
    }

 	 
	/* 
	 * Function for running fastQC on input samples
 	 */
	
    static def fastqc(read, cpus, debug="no") {

    """
        if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi
  		\$print fastqc ${read} 
    """
	}




}