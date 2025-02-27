/*
* NanoPolish
*/

params.LABEL = ""
params.LABELPREP = ""
params.LABELINF = ""
params.EXTRAPARS = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.CONTAINER = 	"quay.io/biocontainers/m6anet:2.1.0--pyhdfd78af_0"

include { unzipCmd } from '../global_functions.nf'

/*
*/
process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
		m6anet -v       
    """
}

process dataprep {
//    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    container params.CONTAINER
    label (params.LABEL)
    label (params.LABELPREP)
    tag "${id}" 
 		
    input:
    tuple val(id), path(eventalign)
    
    output:
    tuple val(id), path("${id}")
    
    script:

    def unzip = unzipCmd(eventalign)
    def cmd = unzip[1]
    def name = unzip[0]
    def clean = unzip[2]

	"""
	${cmd} 
    m6anet dataprep --eventalign ${name} \
                --out_dir ${id} --n_processes ${task.cpus}
    ${clean}
    """



}




/*
*/

process inference {
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    container params.CONTAINER
    label (params.LABEL)
    label (params.LABELINF)
    
    tag "${sampleID}" 
 		
    input:
    tuple val(sampleID), path(input)
    
    output:
    tuple val(sampleID), path("${sampleID}_inf")
    
    script:
	"""
    m6anet inference --input_dir ./${input} --out_dir ./${sampleID}_inf ${params.EXTRAPARS}
    for i in ${sampleID}_inf/*.csv; do gzip \$i; done
    """



}

/*
*/

workflow INFERENCE_RNA004 {

    take: 
    input
    
    main:
    m6data = dataprep(input)
    
	out = inference(m6data)

	emit:
	m6data
 }


/*
*/

workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
