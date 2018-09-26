/* 
 * Class for feature counter
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class FeatureCounter {

    /*
     * Properties definition
     */
    
     String input = ''
     String output = ''
     String annotation = ''
     String memory = 0
     String strand = ''
     Integer cpus = 1
     String extrapars = ''

    /* 
     *  count feature with htseq-count
     */ 
     
    def public htseqCount() {

    """
        htseq-count ${this.extrapars} -s ${this.strand} -f bam ${this.input} ${this.annotation} ${this.extrapars} > ${this.output}
    """
    }

 

}
