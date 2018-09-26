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
     String chr_sizes = ''
     String memory = ''
     Integer cpus = 1
     String extrapars = ''


    def public doPeakcall(peakcaller) { 
     switch (peakcaller) {
        case "macs2":
            this.peakCallWithMacs2()
            break
       case "epic":
            this.peakCallWithEpic()
            break
        default:
            break
        }   
    }
    
    
    /* 
     *  Peak call with MACS
     */ 
    def private peakCallWithMacs2() {

    """
    export PYTHONPATH=""
    macs2 callpeak \
        -t ${this.sample} \
        -f BAM \
        -n ${this.id} \
        --extsize ${this.frag_len} \
        ${extrapars} \
        -g ${this.genome_size} 2>${this.id}_log.txt
    """
    }

    /* 
     *  Peak call with EPIC
     */ 
    def private peakCallWithEpic() {

    """
        epic --treatment ${this.sample} \
             --control ${this.control} \
             --chromsizes ${this.chr_sizes} \
             --effective-genome-fraction ${this.genome_size} \
             -o ${this.id}.out \
             -b ${this.id}.bed \
             --fragment-size ${this.frag_len} \
             ${extrapars}
    """
    }

}
