/* 
 * Repository for biological functions
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class Trimming {

	/* 
	 * Function for trimming read pairs with Skewer
 	 */
	
    static def trimPairsWithSkewer( pair_id, reads, minread, cpus, extrapars="", debug="no") { 
	    """
		if [ `echo ${debug} == "debug"` ]; then print="echo "; else print=""; fi
		\$print	skewer ${extrapars} -t ${cpus} -l ${minread} -n -u -o ${pair_id} -z ${reads}
        """
     }

}