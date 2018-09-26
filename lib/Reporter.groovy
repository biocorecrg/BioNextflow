/* 
 * Class for making reports, file descriptors, etc  
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class Reporter {

    /*
     * Properties definition
     */
    
     String email = ''
     String id = ''
     String title = ''
     String subtitle = ''
     String strand = ''
     String extension = ''
     String type = ''
     String file_name = ''
     String PI = ''
     String user = ''
     String config_file = ''
     String application = ''
     

    def public test() {
        output = this.dump()
    """
        echo '${output}'
    """
    }

   /* 
    * Method for writing files needed for an HUB @ UCSC Genome Browser
    */
    
    def public makeGenomeUcscHub() {
    
    """
        echo "genome ${this.id}
trackDb ${this.id}/trackDb.txt" > genomes.txt;
        echo "hub ${this.id}
shortLabel ${this.title}
longLabel ${this.subtitle}
genomesFile genomes.txt
email ${this.email}" > hub.txt 
    """
    }

   /* 
    * Method for writing track DB file needed for HUB @ UCSC Genome Browser
    */
    def public makeBigWigTrackDB() {

        """
         echo "track ${this.id} 
type ${this.type}  
compositeTrack on 
shortLabel RNAseq profiles (${this.title})
longLabel RNAseq profiles (${this.subtitle})
visibility full
color 0,0,255
autoScale on
maxHeightPixels 128:60:11
" >>  trackDb.txt;
    for i in *.${this.extension}; \\
        do NAME=`echo \$i | sed s/\\.${this.extension}//g`; \\
        echo track \$NAME; \\
        echo bigDataUrl \$i; \\
        echo shortLabel \$NAME; \\
        echo longLabel \$NAME; \\
        echo "type ${this.type}"; \\
        echo "parent ${this.id}"; \\
    echo ""; \\
         done >> trackDb.txt;
    """
    }
    
    
    /*
     * write multiQC config file (config.yaml)
     */

    def public assembleConfigForMultiQC() {
        def myString = "title: \"${this.title}\"\nsubtitle: \"${this.subtitle}\"\nintro_text: False\n\nreport_header_info:\n"
        if (this.PI) { myString +=  "    - PI: ${this.PI}\n" } 
        if (this.user) { myString +=  "    - User: ${this.user}\n" } 
        myString +=  "    - Date: `date`\n    - Contact E-mail: \'${this.email}\'\n    - Application Type: \'${this.application}\'\n    - Reference Genome: \'${this.id}\'\n"
        return (myString)
}
    
   /* 
    * Method for writing multiQC report in a super custom way
    */
    def public makeMultiQCreport() {
        def myString = this.assembleConfigForMultiQC()
        """
            cat > config.yaml << EOL
${myString}
EOL
            cat ${this.config_file} >> config.yaml;   
            multiqc -c config.yaml .
        """
    }

}