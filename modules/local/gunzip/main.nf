// unzipping
process EVENTUAL_UNZIP {
    label ("process_medium")
    tag "${meta.id}"

    container 'quay.io/biocontainers/gzip:1.11'

    input:
    tuple val(meta), path (file, name: "my-dir/*")

    output:
    tuple val(meta), path("*", type:"file")
    
    script:
    if (file.extension == "gz")  {
		"""
 	   	zcat ${file} > ${file.baseName}
 	    """
 	}
    else {
    	"""
    		ln -s ${file} ./${file.name}
    	"""
    }

}
