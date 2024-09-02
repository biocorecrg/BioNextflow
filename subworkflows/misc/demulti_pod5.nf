params.OUTPUT = ""
params.OUTPUTST = ""
params.OUTPUTMODE = "copy"
params.TYPE = "guppy"
params.LABEL = ""

process preparing_demultiplexing_pod5 {

    label (params.LABEL)
    tag "${ idfile }"
		
	input:
	tuple val(idfile), path("demux_*")

	output:
	tuple val(idfile), path("*.files")

	
	script:
	"""
       zcat demux_*|awk 'FNR <= 1' > dem.files
       zcat demux_*| grep -v "read_id"|sort|uniq >> dem.files
	"""
}


process extracting_demultiplexed_pod5_seqtagger {
    label (params.LABEL)
    tag "${ idfile }"

    container = "biocorecrg/pod5:0.2.4"
    
    publishDir(params.OUTPUT, mode:params.OUTPUTMODE, pattern: '*-*')    
   // publishDir(params.OUTPUTST, mode:params.OUTPUTMODE, pattern: 'summaries/*_final_summary.stats', saveAs: { file -> "${file.split('\\/')[-1]}" })    
	
	input:
	tuple val(idfile), path(idlist), file("*")

	output:
	path("${idfile}---*"), type: "dir", emit: dem_fast5
		
	script:
	"""
    pod5 subset *.pod5 -f -t ${task.cpus} --summary ${idlist} --template "${idfile}---barcode-{barcode}.pod5" --columns barcode -M --output ./
	for i in *barcode*.pod5; do mkdir \$(basename "\$i" ".pod5"); cp \$i  \$(basename "\$i" ".pod5"); done
	"""
} 

process extracting_demultiplexed_pod5_dorado {
    label (params.LABEL)
    tag "${ idfile }"

    container = "biocorecrg/pod5:0.2.4"
    
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


workflow DEMULTI_POD5 {

    take: 
    input_stats
    input_pod5    
    
    main:
       prep_demux = preparing_demultiplexing_pod5(input_stats)
       pod5s = input_pod5.transpose().groupTuple()
	   input_data = prep_demux.transpose().combine(pod5s,  by: 0)
	   
       switch(params.TYPE) {                      
          case "dorado":
   			    extracting_demultiplexed_pod5_dorado(input_data)            
          break;
          case "seqtagger":
   			    extracting_demultiplexed_pod5_seqtagger(input_data)            
          break;  
    	}
 
 }
 
 workflow DEMULTI_POD5_FILTER {

    take: 
    input_stats
    input_pod5 
    barcodes   
    
    main:
       prep_demux = preparing_demultiplexing_pod5(input_stats)
       pod5s = input_pod5.transpose().groupTuple()
	   filt_prep_demux = filterDemuxBarcodes(prep_demux.combine(barcodes))
	   input_data = filt_prep_demux.transpose().combine(pod5s,  by: 0)
       
       switch(params.TYPE) {                      

          case "dorado":   
   			    extracting_demultiplexed_pod5_dorado(input_data)            
          break;
          case "seqtagger":
   			    extracting_demultiplexed_pod5_seqtagger(input_data)            
          break;  
    	}
 
 }

process filterDemuxBarcodes {

    tag "${ idfile }"
	
	input:
	tuple val(idfile), path(dem_files), path(barcodes)

	output:
	tuple val(idfile), path("new.dem.files")

	script:
	"""
		awk -F "---" '{print \$2}' ${barcodes} > selected_barcodes.txt
		awk 'FNR <= 1' ${dem_files} > new.dem.files
		awk '{OFS="\t"; split(\$3,a,"_"); print \$1,\$2,a[2]}' ${dem_files} | grep -f selected_barcodes.txt -w >> new.dem.files
	"""
}


