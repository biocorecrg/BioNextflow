/* 
 * Class for NGS mappers 
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 */
 
 class NGSaligner {

    /*
     * Properties definition
     */
     String reference_file = ''
     String annotation_file = ''
     String id = ''
     String index = ''
     String dbtype = ''
     String lib_type = ''
     String reads = ''
     Integer read_size = 0
     String read1 = ''
     String read2 = ''
     String output = ''
     String result = ''
     Integer cpus = 1
     Integer memory = 2
     Integer mismatches = 0
     String extrapars = ''
    
    /*
     * Methods definition
     */
     
    def public test() {
        output = this.dump()
    """
        echo '${output}'
    """
    }

    /*
     * unique interface for alignment
     */
     
    def public doAlignment(String aligner) { 
     switch (aligner) {
        case "STAR":
            this.alignWithStar()
            break
        case "bowtie2":
            this.alignWithBowtie2()
            break
        case "bowtie":
            this.alignWithBowtie()
            break
        case "shortStack":
            this.alignWithShortStack()
            break
        case "blast":
            this.alignWithBlast()
            break
        default:
            break
        }   
    }

    /*
     * unique interface for indexing genomes / transcriptomes
     */ 

    def public doIndexing(String aligner) { 
     switch (aligner) {
        case "STAR":
            this.indexWithSTAR()
            break
         case "bowtie2":
            this.indexWithBowtie2()
            break
         case "bowtie":
            this.indexWithBowtie()
            break
        case "blast":
            this.indexWithBlast()
            break   
        default:
            break
        }   
    }

    /*
     * unique interface for extracting information from indexes
     */ 
    def public genomeStatsFromIndex(String aligner) { 
     switch (aligner) {
         case "bowtie2":
            this.genomeStatsFromBowtie2Index()
            break
         case "bowtie":
            this.genomeStatsFromBowtieIndex()
            break      
        default:
            break
        }   
    }

/*******************************************************************
* INDEXING
********************************************************************/
   /*
    * Indexing a genome with STAR mapper. It reads both gzipped and plain fasta
    */ 
    
    def private indexWithSTAR() {
        """
        annotation_cmd=""
        if [[ `echo ${this.annotation_file}` != "" ]]; then
	    annotation_cmd="--sjdbOverhang ${this.read_size} --sjdbGTFfile ${this.annotation_file}"
        fi
        mkdir ${this.index}
        if [ `echo ${this.reference_file} | grep ".gz"` ]; then 
            zcat ${this.reference_file} > `basename ${this.reference_file} .gz`
            STAR --runMode genomeGenerate --genomeDir ${this.index} --runThreadN ${this.cpus} \
            --genomeFastaFiles `basename ${this.reference_file} .gz` \$annotation_cmd \
            --outFileNamePrefix ${this.index} \
            ${this.extrapars};
            rm `basename ${this.reference_file} .gz`
        else 
            STAR --runMode genomeGenerate --genomeDir ${this.index} --runThreadN ${this.cpus} \
            --genomeFastaFiles ${this.reference_file} \$annotation_cmd \
            --outFileNamePrefix ${this.index} \
            ${this.extrapars}
        fi
        """
    }


    /*
     * Mapping SE and PE reads with STAR. Reads can be both gzipped and plain fastq
     */ 
     
    def private alignWithStar(){
        """
        if [ `echo "${this.reads}"| cut -f 1 -d " " | grep ".gz"` ]; then gzipped=" --readFilesCommand zcat "; else gzipped=""; fi
            STAR --genomeDir ${this.index} \
                 --readFilesIn ${this.reads} \
                  \$gzipped \
                  --outSAMunmapped None \
                  --outSAMtype BAM SortedByCoordinate \
                  --runThreadN ${this.cpus} \
                  --outFileNamePrefix ${this.id} \
                  ${this.extrapars}

                  mkdir ${this.output}
                  mv ${this.id}Aligned* ${this.output}/.
                  mv ${this.id}SJ* ${this.output}/.
                  mv ${this.id}Log* ${this.output}/.
                  if test -f "${this.id}ReadsPerGene*"; then
                      mv ${this.id}ReadsPerGene* ${this.output}/.
                  fi      
        """
    }

    /*
     * Indexing genomes with Bowtie2. Sequence can be both gzipped and plain fasta
     */

    def private indexWithBowtie2() { 
 
        """
        if [ `echo ${this.reference_file} | grep ".gz"` ]; then 
            zcat ${this.reference_file} > `basename ${this.reference_file} .gz`
            bowtie2-build --threads ${this.cpus} `basename ${this.reference_file} .gz` ${this.index} ${this.extrapars}
            rm `basename ${this.reference_file} .gz`
        else bowtie2-build --threads ${this.cpus} ${this.reference_file} ${this.index} ${this.extrapars}
        fi
        """
    }
    
    /*
     * Indexing genome with Bowtie. Sequence can be both gzipped and plain fasta
     */ 
     
    def private indexWithBowtie() { 

        """ 
        if [ `echo ${this.reference_file} | grep ".gz"` ]; then 
            zcat ${this.reference_file} > ${this.index}.fa          
            bowtie-build --threads ${this.cpus} ${this.index}.fa ${this.index} ${this.extrapars}
        else ln -s ${this.reference_file} ${this.index}.fa 
            bowtie-build --threads ${this.cpus} ${this.index}.fa ${this.index} ${this.extrapars}
        fi
        """
    }
    
    /*
     * Indexing genome / proteome with Blast. Sequence can be both gzipped and plain fasta
     */ 
     
    def private indexWithBlast() { 
 
        """ 
        if [ `echo ${this.reference_file} | grep ".gz"` ]; then 
            zcat ${this.reference_file} > ${this.index}         
        else ln -s ${this.reference_file} ${this.index} 
        fi
        makeblastdb -in ${this.index} -dbtype ${this.dbtype} ${this.extrapars}
        """
    }
 
/*******************************************************************
* MAPPING
********************************************************************/

    /*
     * Mapping SE and PE reads with Bowtie2. Reads can be both gzipped and plain fastq
     */ 
    def private alignWithBowtie2() { 
        """
        if [[ "${this.read1}" == "" &&  ${this.reads} != "" ]]; then 
            bowtie2 --non-deterministic -x ${this.index} -U ${this.reads} -p ${this.cpus} ${this.extrapars} | samtools view -Sb -@ ${this.cpus} - > ${this.output}
        else 
            bowtie2 -x ${this.index} -1 ${this.read1} -2 ${this.read2} -p ${this.cpus} ${this.extrapars} | samtools view -Sb -@ ${this.cpus} - > ${this.output}         
        fi  
        """
    }

    /*
     * Mapping SE and PE reads with Bowtie. Reads can be both gzipped and plain fastq
     */ 

    def private alignWithBowtie() { 
        """
        if [[ "${this.read1}" == "" &&  ${this.reads} != "" ]]; then 
            bowtie --non-deterministic -x ${this.index} -U ${this.reads} -p ${this.cpus} ${this.extrapars} | samtools view -Sb -@ ${this.cpus} - > ${this.output}
        else 
            bowtie -x ${this.index} -1 ${this.read1} -2 ${this.read2} -p ${this.cpus} ${this.extrapars} | samtools view -Sb -@ ${this.cpus} - > ${this.output}          
        fi  
        """
    }

    /*
     * Mapping SE reads with ShortStack. This tool is wraps bowtie1.
     */ 

    def private alignWithShortStack() { 
        """
            ShortStack ${this.extrapars} --mismatches ${this.mismatches} --nohp --readfile ${this.reads} --outdir ${this.output} --genomefile ${this.index} --bowtie_cores ${this.cpus}
            cp ${this.output}/ErrorLogs.txt ${this.output}/${id}.txt
        """
    }
    
    /*
     * Aligning sequences with Blast. Sequences can be both gzipped orr plain fasta file
     */ 

    def private alignWithBlast() { 
        """
            blastn -out ${this.output} -db ${this.index} -query ${this.reads} -num_threads ${this.cpus} ${this.extrapars}
        """
    }       
    
/*******************************************************************
* OTHER
********************************************************************/

    /*
     * Get genome STATS from Bowtie2 index
     */ 

    def public genomeStatsFromBowtie2Index() { 
        """
           bowtie2-inspect --summary ${this.index} | awk -F"\t" '{if (\$1~"Sequence") {split(\$2, a, " "); print a[1]"\t"\$3}}' > ${this.output}
        """
    }
    
    /*
     * Get genome STATS from Bowtie index
     */ 

    def private genomeStatsFromBowtieIndex() { 
        """
           bowtie-inspect --summary ${this.index} | awk -F"\t" '{if (\$1~"Sequence") {split(\$2, a, " "); print a[1]"\t"\$3}}' > ${this.output}
        """
    }

}
