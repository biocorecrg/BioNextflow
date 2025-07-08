params.OUTPUT = ""
params.OUTPUTST = ""
params.OUTPUTMODE = "copy"
params.TYPE = "guppy"
params.LABEL = ""
params.CONTAINER = "biocorecrg/pod5:0.2.4"
params.MAX_PROC = 100

process preparing_demultiplexing_pod5 {

    label (params.LABEL)
    tag "${ idfile }"
    container = params.CONTAINER
		
	input:
	tuple val(idfile), path("demux_*")

	output:
	tuple val(idfile), path("*.files")

   script:
   switch(params.TYPE) {                      
	
   	case "seqtagger":   
		"""
       zcat demux_*|awk 'FNR <= 1' > dem.files
       zcat demux_*| grep -v "read_id"|sort|uniq| awk -F "\t" '{OFS="\t"; print \$1,\$2,"bc_"\$3,\$4,\$5 }' >> dem.files
		"""
    break;
    case "dorado":
	"""
       zcat demux_*|awk 'FNR <= 1' > dem.files
       zcat demux_*| grep -v "read_id"|sort|uniq >> dem.files
	"""
    break;  
 }
}



process extracting_demultiplexed_pod5 {
    label (params.LABEL)
    tag "${ idfile }"

    container = params.CONTAINER
    
    publishDir(params.OUTPUT, mode:params.OUTPUTMODE, pattern: '*-*')    
   // publishDir(params.OUTPUTST, mode:params.OUTPUTMODE, pattern: 'summaries/*_final_summary.stats', saveAs: { file -> "${file.split('\\/')[-1]}" })    
	
	input:
	tuple val(idfile), path(idlist), file("*")

	output:
	path("${idfile}---*"), type: "dir", emit: dem_fast5
		
	script:
	"""
    pod5 subset *.pod5 -f -t ${task.cpus} --summary ${idlist} --template "${idfile}---{barcode}.pod5" --columns barcode -M --output ./temp_out
	for i in temp_out/*.pod5; do mkdir \$(basename "\$i" ".pod5"); cp \$i  \$(basename "\$i" ".pod5"); done
	"""
} 

/**
**/

process extracting_demultiplexed_pod5_dorado {
    label (params.LABEL)
    tag "${ idfile }"

    container = params.CONTAINER
    
    publishDir(params.OUTPUT, mode:params.OUTPUTMODE, pattern: '*---*')    
	
	input:
	tuple val(idfile), path(idlist), file("*")

	output:
	path("${idfile}---*"), type: "dir", emit: dem_fast5
		
	script:
	"""
    pod5 subset *.pod5 -f -t ${task.cpus} --summary ${idlist} --template "${idfile}---{barcode}.pod5" --columns barcode -M --output ./
	for i in ${idfile}---*.pod5; do mkdir -p \$(basename "\$i" ".pod5"); mv \$i  \$(basename "\$i" ".pod5"); done
	"""
} 


process split_readids {
	container "biocorecrg/pod5:0.2.4"
    label (params.LABEL)
    tag { sampleID }
    container = params.CONTAINER


    input:
    tuple val(sampleID), path(pod5), val(num)
    
    output:
    tuple val(sampleID), path("pieces.*")

    script:
	"""
		pod5 view -t ${task.cpus} -I -H ${pod5} > read.list
		split -d -l ${num} read.list pieces.
	"""
}

process split_pod5 {
    container = params.CONTAINER
    label (params.LABEL)
    tag { sampleID }
    publishDir(params.OUTPUT, mode:'copy')
    maxForks params.MAX_PROC 

    input:
    tuple val(sampleID), path(pieces), path(pod5)
    
    output:
    tuple val(sampleID), path("*.pod5")

    script:
    def outfile = "${sampleID}_${pieces}.pod5"
	"""
		pod5 filter -t ${task.cpus}  -M -D -i ${pieces} -o ${outfile} ${pod5}
	"""
}

workflow SPLIT_POD5 {

    take: 
    input_pod5
    read_num    
    
    main:
 		read_ids = split_readids(input_pod5.combine(read_num))
		pod5_to_be_split = read_ids.transpose().combine(input_pod5, by: 0)
		split_pod5(pod5_to_be_split)
	
 }
 

workflow DEMULTI_POD5 {

    take: 
    input_stats
    input_pod5    
    
    main:
       prep_demux = preparing_demultiplexing_pod5(input_stats)
       pod5s = input_pod5.transpose().groupTuple()       
	   input_data = prep_demux.transpose().combine(pod5s,  by: 0)
	   extracting_demultiplexed_pod5(input_data)
  
 }
 
 workflow DEMULTI_POD5_FILTER {

    take: 
    input_stats
    input_pod5 
    barcodes   
    
    main:
       prep_demux = preparing_demultiplexing_pod5(input_stats)
       pod5s = input_pod5.transpose().groupTuple()       
	   filt_prep_demux = filxterDemuxBarcodes(prep_demux.combine(barcodes, by: 0))
	   input_data = filt_prep_demux.transpose().combine(pod5s,  by: 0)
   	   extracting_demultiplexed_pod5(input_data)            
       
 
 }

process filterDemuxBarcodes {
    container = params.CONTAINER
	
    tag "${ idfile }"
	
	input:
	tuple val(idfile), path(dem_files), path(barcodes)

	output:
	tuple val(idfile), path("new.dem.files")

	script:
    switch(params.TYPE) {                      

    case "seqtagger":   
		"""
		awk -F "---" '{print \$2}' ${barcodes} > selected_barcodes.txt
		awk 'FNR <= 1' ${dem_files} > new.dem.files
		awk '{OFS="\t"; split(\$3, a, "_"); print \$1, \$2, "bc_"(length(a) > 1 ? a[2] : \$3)}' ${dem_files} | grep -f selected_barcodes.txt -w >> new.dem.files
		"""
    break;
    case "dorado":
	"""
		awk -F "---" '{print \$2}' ${barcodes} > selected_barcodes.txt
		awk 'FNR <= 1' ${dem_files} > new.dem.files
		awk '{OFS="\t"; split(\$3, a, "_"); print \$1, \$2, (length(a) > 1 ? a[2] : \$3)}' ${dem_files} | grep -f selected_barcodes.txt -w >> new.dem.files
	"""
    break;  
    }
 }



