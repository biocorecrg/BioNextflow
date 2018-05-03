/* 
 * Repository for biological functions
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class QC {

	/* 
	 * Function for qualimap QC

	
    static def qualimapRNAseq ( bamfile, annotation_file, outfolder="QUALIMAP", strand="strand-specific-reverse", memory="2G", read_type="SINGLE", extrapars="", debug="no") { 

    """
       if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi
       if [ `echo ${single} == "SINGLE"` ]; then pe_mode="-pe"; else if [`echo ${single} == "PAIRED"`] pe_mode=""; fi; fi
      
	   \$print unset DISPLAY
       \$print mkdir tmp
       \$print export JAVA_OPTS="-Djava.awt.headless=true -Xmx${memory} -Djava.io.tmpdir=\$PWD/tmp"    
       \$print qualimap rnaseq \${pe_mode} --java-mem-size=${memory} -bam ${bamfile} -gtf ${annotation_file} -outdir ${outfolder} -p ${strand} ${extrapars}
       \$print rm -fr tmp
    """
    }
	*/
 	 
	/* 
	 * Function for running fastQC on input samples
 	 */
	
    static def fastqc(read, cpus="1", debug="no") {

    """
        if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi
  		\$print fastqc -t ${cpus} ${read} 
    """
	}

	/* 
	 * Function for running fastQC on input samples
	 | 
 	 */
	
    static def getReadSize(read, debug="no") {

    """
        if [ `echo ${read} | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi

        if [ `echo ${debug} == "debug"` ]; then
			echo "complex to escape!" 
        else
		\$cat ${read} | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", sum/100; exit} }'
		fi
    """
	}



}