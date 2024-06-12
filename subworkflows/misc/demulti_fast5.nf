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

process preparing_demultiplexing_fast5_deeplexicon {

    tag "${ idfile }"
		
	input:
	tuple val(idfile), path("demux_*")

	output:
	tuple val(idfile), path("*.list")

	
	script:
	"""
	cat demux_* | grep -v ReadID >> dem.files
	awk '{print \$2 > \$3".list" }' dem.files
	"""
}

process preparing_demultiplexing_fast5_seqtagger {

    label (params.LABEL)
    tag "${ idfile }"
		
	input:
	tuple val(idfile), path("demux_*")

	output:
	tuple val(idfile), path("*.files")

	
	script:
	"""
	zcat demux_* | sed -e '2,\${/^read_id/d' -e '}' - >> dem.files
	"""
}

process filterDemuxBacodes_seqtagger {
	
	input:
	tuple val(idfile), path(dem_files)
	val barcodes

	output:
	tuple val(idfile), path("filtered_dem.files")

	script:
	"""
	for i in $barcodes; do barcode=\$(cut -d "_" -f 1 <(echo \$i| rev)); awk -v bar=\$barcode 'BEGIN{print "read_id\tadapter_end\tbarcode\tmapQ\tbaseQ"} {if (\$5>50 && \$3==bar){print \$0}}' $dem_files >> demux_file.txt ; done

	sed -e '2,\${/^read_id/d' -e '}' demux_file.txt > filtered_dem.files	
	"""
}

process extracting_demultiplexed_fast5 {
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

process extracting_demultiplexed_fast5_seqtagger {
    label (params.LABEL)
    
	container 'lpryszcz/seqtagger:1.0a3'
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
	mkdir fast5_files_input
	mv *.fast5 ./fast5_files_input

	fast5_split_by_barcode.py -b 50 -f ./fast5_files_input -i $idlist -o ./

	for dir in ./bc* ; do mv \$dir ${idfile}---\$(basename "\$dir") ; done

	mkdir summaries
	for i in */filename_mapping.txt; do mv \$i summaries/\$(cut -d "/" -f1 <(echo \$i))_final_summary.stats; done
	"""
} 




workflow DEMULTI_FAST5 {

    take: 
    input_stats
    input_fast5    
    
    main:
       switch(params.TYPE) {                      
           case "guppy":   
				input_data = input_stats.join(input_fast5)
                extracting_demultiplexed_fast5_guppy(input_data)
          break;
          case "deeplexicon":
     		    prep_demux = preparing_demultiplexing_fast5_deeplexicon(input_stats)
			    input_data = prep_demux.transpose().combine(input_fast5,  by: 0)
     			extracting_demultiplexed_fast5(input_data)            
          break;  
          case "seqtagger":
      		    prep_demux = preparing_demultiplexing_fast5_seqtagger(input_stats)
			    input_data = prep_demux.transpose().combine(input_fast5,  by: 0)
     			extracting_demultiplexed_fast5_seqtagger(input_data)            
          break;  
    	}
 
 }
 
 workflow DEMULTI_FAST5_FILTER {

    take: 
    input_stats
    input_fast5 
    barcodes   
    
    main:
       switch(params.TYPE) {         
           // not so clear how to filter             
          case "guppy":   
				input_data = input_stats.join(input_fast5)
                extracting_demultiplexed_fast5_guppy(input_data)
          break;
          case "deeplexicon":
     		    prep_demux = preparing_demultiplexing_fast5_deeplexicon(input_stats)
			    filt_prep_demux = filterDemuxBacodes(prep_demux, barcodes)
			    input_data = filt_prep_demux.transpose().combine(input_fast5,  by: 0)
     			extracting_demultiplexed_fast5(input_data)            
          break;  
          case "seqtagger":
      		    prep_demux = preparing_demultiplexing_fast5_seqtagger(input_stats)
			    //filt_prep_demux = filterDemuxBacodes(prep_demux, barcodes)
				filt_prep_demux = filterDemuxBacodes_seqtagger(prep_demux, barcodes)
			    input_data = filt_prep_demux.transpose().combine(input_fast5,  by: 0)
     			extracting_demultiplexed_fast5_seqtagger(input_data)            
          break;  
    	}
 
 }
 
 // Create a channel for included ids
def filterDemuxBacodes (mylists, mybarcodes) {
	reshaped_barcoded_data = mylists.transpose().map{
		[ "${it[0]}---${it[1].simpleName}", it[1] ] 
	}
		
	filtered_data = reshaped_barcoded_data.combine(mybarcodes,  by: 0).map{
		def id = it[0].split("---")[0]
		[id, it[1]]
	}.groupTuple()
	
	return(filtered_data)
}

