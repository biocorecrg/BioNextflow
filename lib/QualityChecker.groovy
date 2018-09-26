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
    
     String id = ''
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

    /* 
     * Function to calculate fingerprints using deeptools
     * Input is a list of bam files / Output a prefix
     */
    
    def public plotFingerprintWithDeepTools() {

    """
        plotFingerprint \
         -b ${this.input} \
        --smartLabels \
        --minMappingQuality 30 --skipZeros \
        -T "Fingerprints"  \
        --plotFile ${this.output}.png \
        --outRawCounts ${this.output}_rawcounts.tab \
        --numberOfProcessors ${this.cpus} ${this.extrapars} \
    """
    }

    /* 
     * Function to calculate multiBamSummary using deeptools
     * Input is a list of bam files / Output a pnz file
     */
    
    def public multiBamSummaryWithDeepTools() {

    """
        multiBamSummary \
        bins \
         -b ${this.input} \
        --smartLabels \
        --outFileName ${this.output} \
        --numberOfProcessors ${this.cpus} ${this.extrapars} \
    """
    }
    
    /* 
     * Function to plot PCA using deeptools
     * Input is a file obtained by multiBigwigSummary or multiBamSummary
     */
    
    def public plotPCAWithDeepTools() {

    """
        plotPCA \
         -in ${this.input} \
         -T "PCA of read counts" \
        -o ${this.output} \
        ${this.extrapars}
    """
    }

    /* 
     * Function to plot correlation using deeptools
     * Input is a file obtained by multiBigwigSummary or multiBamSummary
     */
    
    def public plotPCorrelationWithDeepTools() {

    """
        plotCorrelation \
         -in ${this.input} \
          --plotTitle "Spearman Correlation of Read Counts" \
          --corMethod spearman --skipZeros \
          --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
          -o ${this.output} \
        ${this.extrapars}
    """
    }

    /* 
     * Function to calc Matrix for gene body using deeptools
     * Input is list of bigwig files
     */
    
    def public computeMatrixForGenesWithDeepTools() {

    """ 
    computeMatrix scale-regions -S ${this.input} \
        -R ${this.annotation_file} \
        --beforeRegionStartLength 3000 \
        --smartLabels \
        --regionBodyLength 5000 \
        --afterRegionStartLength 3000 \
        --numberOfProcessors ${this.cpus} \
        --skipZeros -o ${this.output} ${this.extrapars}
    """
    }

    /* 
     * Function to calc Matrix for TSS using deeptools
     * Input is list of bigwig files
     */
    
    def public computeMatrixForTSSWithDeepTools() {

    """ 
    computeMatrix reference-point -S ${this.input} \
        -R ${this.annotation_file} \
        --referencePoint TSS \
        --smartLabels \
        --beforeRegionStartLength 3000 \
        --afterRegionStartLength 3000 \
        --numberOfProcessors ${this.cpus} \
        --skipZeros -o ${this.output} ${this.extrapars}
    """
    }

    /* 
     * Function to plot Heatmap using deeptools
     * Input is list of matrix created by computeMatrix
     */ 
    def public plotHeatmapWithDeepTools() {

    """ 
    plotHeatmap \
        -m ${this.input} \
        -out ${this.output} \
        --colorMap jet \
        --missingDataColor "#FFF6EB" \
        --heatmapHeight 15 \
        ${this.extrapars}
    """
    }

    /* 
     * Function to plot profiles using deeptools
     * Input is list of matrix created by computeMatrix
     */ 
    def public plotProfileWithDeepTools() {

    """ 
    plotProfile -m ${this.input} \
        -out ${this.output} \
        ${this.extrapars}
    """
    }
}
