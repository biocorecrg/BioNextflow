/* 
 * Class for misc functions (to be splitted?) 
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class PeakCaller {

	/*
	 * Properties definition
	 */
	
     String sample = ''
     String control = ''
     String output = ''
     String mode = ''
     String id = ''
     String frag_len = ''
     String genome_size = ''
     String memory = ''
     Integer cpus = 1
     String extrapars = ''


    def public doPeakcall(peakcaller) { 
     switch (peakcaller) {
        case "macs2":
			this.peakCallWithMacs2()
            break
        default:
            break
    	}	
	}
	
	
	/* 
	 *  Sorting bam files with samtools
 	 */	
    def private peakCallWithMacs2() {

    """
    export PYTHONPATH=""
    macs2 callpeak \\
        -t ${this.sample} \\
        -f BAM \\
        -n ${this.id} \\
        --fix-bimodal \\
        --extsize ${this.frag_len} \\
        -B \\
        --SPMR \\
        ${extrapars} \\
		-g ${this.genome_size} 2>${this.id}_log.txt
    """
	}

 

}
