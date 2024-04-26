/*
* The accessible subworkflows are:
* GET_VERSION that emits the version of bwa and samtools as stdout*
* The parameters are: 
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for bwa
*	OUTPUT for storing the final sub-workflow output 
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "quay.io/biocontainers/snpsplit:0.6.0--hdfd78af_0"

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		SNPsplit --version| grep Version
    """
}


process index_dual {
    label (params.LABEL)
    tag { "${reference}" }
    container params.CONTAINER

    input:
    tuple path(reference, stageAs: 'genome/*'), path(vcf), val(genomeA),  val(genomeB) 

    output:
    path("genome.masked.fa.gz"), emit: genome
    path("*.based_on*.txt.gz"), emit: snps
    
    script:
    """
    if [ `echo "${reference}" | grep ".gz"` ]; then zcat genome/*.gz >> genome/genome.fa; fi
       SNPsplit_genome_preparation ${params.EXTRAPARS} --vcf_file ${vcf} --reference_genome ./genome/ --dual_hybrid --strain ${genomeA} --strain2 ${genomeB}
	   cat *_N-masked/*.fa | gzip -c >> genome.masked.fa.gz 
	   gzip *.based_on*.txt
    """
}

process split {
    label (params.LABEL)
    tag { pair_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(pair_id), path(reads)
    path(indexes)

    output:
    tuple val(pair_id), path("${pair_id}.bam") 
    
	script:
    def indexname = indexes[0].baseName - ~/\.\w+$/

    """
    biscuit align -@ ${task.cpus} ${params.EXTRAPARS} ${indexname} ${reads} > ${pair_id}.bam
    """
}





workflow SPLIT {
    take: 
    input
    indexes
    
    main:
		out = split(input, indexes)
    emit:
    	out
}

workflow INDEX {
    take: 
    reference
    vcf
    genomeA
    genomeB
    
    main:
		out = index_dual(reference.combine(vcf).combine(genomeA).combine(genomeB))

    emit:
    	genome = out.genome
    	snps = out.genome
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
