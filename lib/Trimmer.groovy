/* 
 * Class for trimming - filtering tools
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class Trimmer {

    /*
     * Properties definition
     */
    
     String reads = ''
     String read = ''
     String read_A = ''
     String read_B = ''
     String id = ''
     Integer min_read_size = 0
     String memory = 0
     Integer cpus = 0
     String extrapars = ''

    /* 
     * Function for trimming read pairs with Skewer
     */
    
    def public trimWithSkewer() { 
        """
            skewer ${this.extrapars} -t ${this.cpus} -l ${this.min_read_size} -n -u -o ${this.id} -z ${this.reads}
        """
     }

}
