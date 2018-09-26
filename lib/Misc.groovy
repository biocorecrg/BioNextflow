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
    
     String id = ''
     String input = ''
     String output = ''
     String mode = ''
     Integer number = 0
     Integer read_size = 0
     String memory = 0
     String genome_size = ''
     String chr_size_file = ''
     String java_path = ''
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
     * Indexing genome with samtools faidx (compressed sequences too)
     */
    
    def public st_indexGenome() { 
    """
        if [ `echo ${this.input} | grep ".gz"` ]; then 
            zcat ${this.input} > `basename ${this.input} .gz`
            samtools faidx `basename ${this.input} .gz`
        else samtools faidx ${this.input}
        fi
    """
    }
    
    /* 
     * Function for getting transcript from GTF 
     */
    static def getTranscriptsFromGTF(genome_file, annotation, output="transcript.fa", debug="no") {

    """
        if [ `echo ${genome_file} | grep ".gz"` ]; then 
            zcat ${genome_file} > `basename ${genome_file} .gz`
            gffread -g `basename ${genome_file} .gz` -w ${output} ${annotation}
            rm `basename ${genome_file} .gz`
        else gffread -g ${genome_file} -w ${output} ${annotation}
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

    /* 
     * Function for making coverage profiles. This is better than the makeAlnProfiles
     * WARNING there is no distinction between strands and it will be replaced in near future with something better
     */

    def public makeAlnProfilesWithDeepTools() {
    """
        bamCoverage --bam ${this.input} -o ${this.output} \
        --binSize 10 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize ${this.genome_size} \
        --extendReads  ${this.read_size}  \
        ${this.extrapars} \
        -p ${this.cpus} 
    """
    }

    /* 
     * Function for estimating Effective Genome Size by using EPIC's script epic-effective
     * You need both EPIC and Jellyfish
     */

    def public estimateEffectiveGenomeSizeWithEpic() {
        """
        if [ `echo ${this.input} | grep ".gz"` ]; then 
            zcat ${this.input} > `basename ${this.input} .gz`
            epic-effective -t ./ --read-length=${this.read_size} -n ${this.cpus} `basename ${this.input} .gz` 2>/dev/null > ${this.output}
            rm `basename ${this.input} .gz`
        else
            epic-effective -t ./ --read-length=${this.read_size} -n ${this.cpus} basename ${this.input} 2>/dev/null > ${this.output}
        fi
        """
    }


    /* 
     * Mark read duplicate using Picard
     * Reads must be sorted!
     */

    def public markDuplicateWithPicard() {
        """
            java -jar ${this.java_path} MarkDuplicates \
            INPUT=${this.input} \
            OUTPUT=${this.output} \
            ASSUME_SORTED=true \
            REMOVE_DUPLICATES=true \
            METRICS_FILE=${this.id}.metrics.txt \
            VALIDATION_STRINGENCY=LENIENT \
            ${this.extrapars}
        """
    }
    
}
