process RIBODETECTOR {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/ribodetector:0.3.1--pyhdfd78af_0' :
        'quay.io/biocontainers/ribodetector:0.3.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_rna_perc.txt")           , emit: output_perc
    tuple val(meta), path("rna*.fq")                  , emit: rrna_fastqs
    path  "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    
    if (meta.single_end) {
        """
   		tot=`awk '{num++}END{print num/4}' ${reads}`
        length=`awk '{num++}{if(num%4==2) {seq++; len+=length(\$0)}}END{print int(len/seq)}' ${reads}`

		ribodetector_cpu  -t ${task.cpus} \
		-l \$length \
		-o /dev/null /dev/null -e rrna -r rna1.fq \
		 -i ${reads} 

		awk -v id=${prefix} -v tot=\$tot '{num++}END{print id" "num/4/tot*100}' rna1.fq > ${prefix}_rna_perc.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ribodetector: \$( ribodetector --version 2>&1 | sed 's/ribodetector //g' )
        END_VERSIONS


        """
        
    } else {
     
		"""
			tot=`awk '{num++}END{print num/4}' ${reads[0]}`
			length=`awk '{num++}{if(num%4==2) {seq++; len+=length(\$0)}}END{print int(len/seq)}' ${reads[0]}`
	
			ribodetector_cpu  -t ${task.cpus} \
			-l \$length \
			-o /dev/null /dev/null -e rrna -r rna1.fq rna2.fq \
			 -i ${reads} 
	
			awk -v id=${prefix} -v tot=\$tot '{num++}END{print id" "num/4/tot*100}' rna1.fq > ${prefix}_rna_perc.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ribodetector: \$( ribodetector --version 2>&1 | sed 's/ribodetector //g' )
        END_VERSIONS

	
         """
    
    }

}
