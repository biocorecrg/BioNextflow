params.OUTPUT = ""
params.OUTPUTST = ""
params.OUTPUTMODE = "copy"
params.TYPE = "guppy"
params.LABEL = ""

process extracting_demultiplexed_fast5_guppy {

    tag "${ idfile }"
    label (params.LABEL)

    publishDir(params.OUTPUT, mode:params.OUTPUTMODE, pattern: '*-*')   
    publishDir(params.OUTPUTST, mode:params.OUTPUTMODE, pattern: 'summaries/*_final_summary.stats', saveAs: { file -> "${file.split('\\/')[-1]}" })    

    container "quay.io/biocontainers/ont-fast5-api:4.0.0--pyhdfd78af_0"
             
	input:
	tuple val(idfile), path("summaries_*"), file("*")
    
	output:
	path("${idfile}-*"), type: "dir", emit: dem_fast5
	path("summaries/*_final_summary.stats"), emit: dem_summaries

    script:
    """
      if [ -f "summaries_" ]; then
	  ln -s summaries_ final_summary.stats
	  else 
		  head -n 1 summaries_1 > final_summary.stats
	      for i in summaries_*; do cat \$i | awk -F"\t" -v id=${idfile} '
	      NR==1 {
    		for (i=1; i<=NF; i++) {
        	  f[\$i] = i
   		    }
          }	      
	      {OFS="\t"; \$(f["barcode_arrangement"]) = id"---"\$(f["barcode_arrangement"]); if (\$1!= "filename") { print \$0} }'  >> final_summary.stats; done
	  fi

		demux_fast5 -c vbz -t ${task.cpus} --input ./ --save_path ./ --summary_file final_summary.stats 

	    mkdir summaries
	    for i in */filename_mapping.txt; do awk 'BEGIN{print "filename\tread_id"}{print \$2"\t"\$1}' \$i > `echo \$i | awk -F"/" '{print "summaries/"\$1"_final_summary.stats"}'`; done

		rm -fr barcode_arrangement
    """
}


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



process extracting_demultiplexed_pod5 {
    label (params.LABEL)
    
	container 'lpryszcz/deeplexicon:1.2.0'
    tag "${ idfile } on ${ idlist }"
    
    publishDir(params.OUTPUT, mode:params.OUTPUTMODE, pattern: '*-*')    
    publishDir(params.OUTPUTST, mode:params.OUTPUTMODE, pattern: 'summaries/*_final_summary.stats', saveAs: { file -> "${file.split('\\/')[-1]}" })    

		
	input:
	tuple val(idfile), path(idlist), file("*")

	output:
	path("${idfile}-*"), type: "dir", emit: dem_fast5
	path("summaries/*_final_summary.stats"), emit: dem_summaries
	
	script:
	"""
	mkdir ${idfile}---`basename ${idlist} .list`; fast5_subset --input ./ --save_path ${idfile}---`basename ${idlist} .list`/ --read_id_list ${idlist} --batch_size 4000 -c vbz -t ${task.cpus}
	mkdir summaries
	for i in */filename_mapping.txt; do awk 'BEGIN{print \$2"\t"\$1}' \$i > `echo \$i | awk -F"/" '{print "summaries/"\$1"_final_summary.stats"}'`; done
	rm */filename_mapping.txt;
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
	for i in *barcode*.pod5; do mkdir \$(basename "\$i" ".pod5"); mv \$i  \$(basename "\$i" ".pod5"); done
	"""
} 

process extracting_demultiplexed_pod5_dorado {
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
    pod5 subset *.pod5 -f -t ${task.cpus} --summary ${idlist} --template "${idfile}---{barcode}.pod5" --columns barcode -M --output ./
	for i in ${idfile}---*.pod5; do mkdir \$(basename "\$i" ".pod5"); mv \$i  \$(basename "\$i" ".pod5"); done
	"""
} 


workflow DEMULTI_POD5 {

    take: 
    input_stats
    input_pod5    
    
    main:
       prep_demux = preparing_demultiplexing_pod5(input_stats)
       pod5s = input_pod5.transpose().groupTuple()
	   prep_demux.transpose().view()
	   pod5s.view()
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
	   filt_prep_demux.transpose().view()
	   pod5s.view()
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


