/* 
 * Class for quality control tools 
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class QualityChecker {

	/*
	 * Properties definition
	 */
	
     String input = ''
     String annotation_file = ''
     String output = ''
     String mode = ''
     String strand = ''
     Integer read_number = 0
     String memory = 0
     Integer cpus = 1
     String extrapars = ''
	 

    def public test() {
        output = this.dump()
    """
     echo '${output}'
    """
    }

   /* 
    * Function for qualimap QC for RNAseq
    */
	
    def public qualimapRNAseq() {
    
    """
	   pe_mode="";
       if [[ "${this.mode}" = "pe" ]]; then pe_mode="-pe"; fi
       
	   unset DISPLAY
       mkdir tmp
       export JAVA_OPTS="-Djava.awt.headless=true -Xmx${this.memory} -Djava.io.tmpdir=\$PWD/tmp"    
       qualimap rnaseq \${pe_mode} --java-mem-size=${this.memory} -bam ${this.input} -gtf ${this.annotation_file} -outdir ${this.output} -p ${this.strand} ${this.extrapars}
       rm -fr tmp
    """
    }

   /* 
    * Function for qualimap QC for bam QC
    */
	
    def public qualimapBamQC() {
    
    """
	   pe_mode="";
       if [[ "${this.mode}" = "pe" ]]; then pe_mode="-pe"; fi
       
	   unset DISPLAY
       mkdir tmp
       export JAVA_OPTS="-Djava.awt.headless=true -Xmx${this.memory} -Djava.io.tmpdir=\$PWD/tmp"    
       qualimap bamqc \${pe_mode} --java-mem-size=${this.memory} -bam ${this.input} -gff ${this.annotation_file} -outdir ${this.output} ${this.extrapars}
       rm -fr tmp
    """
    }
 	 
    /* 
     * Function for running fastQC on input samples
     */
	
    def public fastqc() {

    """
  		fastqc -t ${this.cpus} ${this.input} ${this.extrapars}
    """
	}

      /*
       * check ribosomal content using riboPicker
       */

    def public checkRibo() {

    """
        if [ `echo ${this.input} | grep "gz"` ]; then gzipped=" -gz "; else gzipped=""; fi
		check_rRNA_contam.pl \$gzipped -s ${this.read_number} -f ${this.input} ${this.extrapars} > ${this.output}
    """
	}

      /* 
       * Function for getting the read size on fastq files
       */
	
    def public getReadSize() {

    """
        if [ `echo ${this.input} | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
		\$cat ${this.input} | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", sum/100; exit} }'
    """
	}



}
