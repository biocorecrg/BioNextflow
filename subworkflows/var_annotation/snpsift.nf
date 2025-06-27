/*
* snpeff subworkflows 
*
* The parameters are: 
*	LABEL that allows connecting labels specified in nextflow.config with the subworkflows
*	EXTRAPARS only for mapping step for adding custom command line parameters for star
*	OUTPUT for storing the final sub-workflow output 
*	CONTAINER that can be eventually overridden for feeding a custom container from the main.nf file
*/

params.LABEL = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.CONTAINER = "biocorecrg/snpsift:5.2f"


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    snpsift -version
    """
}

process snpsift_ann {
    label (params.LABEL)
    tag "${id}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.vcf.gz') }

    
    input:
    tuple val(id), path(vcf), path(vcf_ref), path(vcf_ref_idx)

    output:
    tuple val(id), path("*.vcf.gz")
    
    script:
    def name_ref = vcf_ref.simpleName
    """
    snpsift annotate ${params.EXTRAPARS} -XX:+PerfDisableSharedMem ${vcf_ref} ${vcf} > ${id}.on.${name_ref}.vcf
    gzip ${id}.on.${name_ref}.vcf
    """
}


process snpsift_makedb {
    label (params.LABEL)
    tag "${id}"
    container params.CONTAINER

    
    input:
    tuple val(id), path(vcf_ref), path(vcf_ref_idx)

    output:
    tuple val(id), path("*.snpsift.vardb")
    
    script:
    // We need to make all the fields in a variable since it does not work without -fields
    """
FIELDS=\$(zcat ${vcf_ref} | awk 'BEGIN { first_data_line = 1 }
  \$0 !~ /^#/ && first_data_line {
    split(\$NF, a, ";")
    for (i = 1; i <= length(a); i++) {
      split(a[i], b, "=")
      keys[i] = b[1]
    }
    for (i = 1; i <= length(keys); i++) {
      printf "%s", keys[i]
      if (i < length(keys)) printf ","
    }
    print ""
    first_data_line = 0
  }
')

    snpsift -Xmx${task.memory.giga}g annmem \
    -create \
    -dbfile ${vcf_ref} -fields "\$FIELDS" 
    """
    
}

process annotate_mem {
    label (params.LABEL)
    tag "${id}"
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy', pattern: '*.vcf.gz') }

    
    input:
    tuple val(id), path(vcf), path(dbs)

    output:
    tuple val(id), path("*.vcf.gz")
    
    script:
    dbslist = dbs.collect { "-dbfile ${it.toString().replace('.snpsift.vardb', '')}" }.join(' ')

    """
    export JAVA_TOOL_OPTIONS="-XX:-UsePerfData"

    snpsift -Xmx${task.memory.giga}g annmem \
    ${dbslist} \
    ${vcf} | gzip -c > ${id}.ann.ondbs.vcf.gz
    """
    
}

workflow ANN {

    take: 
    vcf
    vcf_ref
    vcf_ref_idx
    
    main:
        out = snpsift_ann(vcf.combine(vcf_ref).combine(vcf_ref_idx))
	emit:
    	out
	
}

workflow ANNOTATE_MEM {

    take: 
    vcf
    vcf_refs
    vcf_ref_idxs
    
    main:
		dbs = snpsift_makedb(vcf_refs.combine(vcf_ref_idxs, by:0))
        db_list = dbs.map{it[1]}.toList().map{
        	[ it ]
        }

        annotate_mem(vcf.combine(db_list))
        
	emit:
    	dbs
	
}

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
