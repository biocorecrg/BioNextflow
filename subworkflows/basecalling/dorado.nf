/*
*  This workflow allows to wrap the DORADO program for basecalling on ONT fast5 data
*  when included you can specify the GPU param ON or OFF for using the GPU
*/

params.LABELBC = ""
params.LABELCONV = ""
params.EXTRAPARS = ""
params.EXTRAPARS_BC = ""
params.EXTRAPARS_DEM = ""
params.OUTPUT = ""
params.OUTPUTMODE = "copy"
params.MOP = ""
params.CONTAINER = "docker pull ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22"

params.GPU = ""
params.DUPLEX = ""

def gpu_cmd = ""
def library_export = ""

if (params.GPU == "OFF") {
	gpu_cmd = '-x "cpu"'
}

my_container = params.CONTAINER 

if (params.GPU == "LOCAL") {
        my_container = ""
}

process getVersion {
    container my_container
    label (params.LABELBC)

    output:
	stdout emit: out    
    
    shell:
    """
		dorado --version 2>&1 | head -n1    
    """
}

process baseCall2Fastq {
    tag { idfile }
    label (params.LABELBC)
	if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*.fastq.gz',  mode: params.OUTPUTMODE ) }

    container my_container
             
    input:
    tuple val(idfile), path(fast5), path(models)
    
    output:
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq

    script:
    if (params.DUPLEX == "") {
    
    """
          dorado basecaller ${gpu_cmd} ${params.EXTRAPARS} --emit-fastq ./ > ${idfile}.fastq
          bgzip -@ ${task.cpus} ${idfile}.fastq
    """
    } else {
    """
          dorado duplex ${gpu_cmd} ${params.EXTRAPARS} --emit-fastq ./ > ${idfile}.fastq
          bgzip -@ ${task.cpus} ${idfile}.fastq
    """
    }
}

process baseCall {
    tag { idfile }
    label (params.LABELBC)

    container my_container

    input:
    tuple val(idfile), path(fast5), path(models)
    
    output:
    tuple val(idfile), path("*.bam"), emit: basecalled_bam

    script:

    script:
    if (params.DUPLEX == "") {
    
    """
          dorado basecaller ${gpu_cmd} --models-directory ./ ${params.EXTRAPARS} ./ > ${idfile}.bam
    """
    } else {
    """
          dorado duplex ${gpu_cmd} --models-directory ./ ${params.EXTRAPARS} ./ > ${idfile}.bam
    """
	}
}

process bam2ModFastq {
    tag { idfile }
    label (params.LABELCONV)
	if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*.fastq.gz',  mode: params.OUTPUTMODE ) }

    container params.CONTAINER 
         
    input:
    tuple val(idfile), path(bam)
    
    output:
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq

    script:

    """
          samtools fastq -@ ${task.cpus} -T MN,MM,ML,mv,pt,ts ${bam} > ${idfile}.fastq
          bgzip -@ ${task.cpus} ${idfile}.fastq
    """
}

process bam2Fastq {
    tag { idfile }
    label (params.LABELCONV)
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*.fastq.gz',  mode: params.OUTPUTMODE ) }

    container params.CONTAINER 

    input:
    tuple val(idfile), path(bam)
    
    output:
    tuple val(idfile), path("*.fastq.gz"), emit: basecalled_fastq

    script:

    """
          samtools fastq -@ ${task.cpus} ${bam} > ${idfile}.fastq
          bgzip -@ ${task.cpus} ${idfile}.fastq
    """
}


process demultiPlex {

    tag { idfile }
    label (params.LABELCONV)

    container my_container
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, pattern: '*.fastq.gz',  mode: params.OUTPUTMODE ) }
   
             
    input:
    tuple val(idfile), path(bam)
    
    output:
    tuple val(idfile), path("*.bam"), emit: demulti_bams
    tuple val(idfile), path("barcoding_summary.txt.gz"), emit: bar_summary

    script:

    """
         dorado demux --emit-summary --threads ${task.cpus} --output-dir ./ --no-classify ${bam}
         gzip barcoding_summary.txt
    """
}

process downloadModel {
    // HACK FOR CLEANUP
    scratch true 
 
    tag { idfile }
    label (params.LABELBC)
    container my_container
    
    input:
    tuple val(idfile), path(bam), path(modelfolder)
    
    output:
    path "${modelfolder}/*", type:'dir'


    script:    
    def down_pars = params.EXTRAPARS.split(" ").find { it.contains('@') }
    def down_pars2 = params.EXTRAPARS.trim().tokenize()[0]
  
    script:
    if (params.DUPLEX == "") {
    
    """
    	echo "here"
       if dorado basecaller ${gpu_cmd} ${down_pars2} --max-reads 1 --models-directory \$PWD/${modelfolder} ./ > test.bam; 
        then
        	echo "Automatic model download succeeded"
        else 
        	echo "Trying the manual download...";
	        dorado download --model ${down_pars} --models-directory \$PWD/${modelfolder}
	    fi
    """
    } else {
    """

	touch stderr.txt
	timeout 1m dorado duplex ${gpu_cmd} ${down_pars2} --models-directory \$PWD/${modelfolder} ./ > test.bam 2>stderr.txt || ( [[ \$? -ne 0 ]] &&  echo "Timeout reached" )

	# Check if the dorado process succeeded or failed
	if grep -q "Starting Stereo Duplex pipeline" stderr.txt; then
    	    echo "Automatic model download succeeded"
	else
        echo "Trying the manual download..."
        dorado download --model ${down_pars} --models-directory \$PWD/${modelfolder}
	fi

	"""
    }
 }


 workflow BASECALL_DEMULTI {
    take: 
    input_fast5
    model_folder
    
    main:
     	model_folders = downloadModel(input_fast5.first().combine(model_folder))
        models = model_folders.collect().map{ [ it ] }
    	bam = baseCall(input_fast5.combine(models))
    	dem_res = demultiPlex(bam)
    	demulti_bams = dem_res.demulti_bams.transpose().map{
            def bam_name = it[1].getSimpleName()
            def barcode = "${bam_name}".split("_").last()
    	    def new_id = "${it[0]}.${barcode}"
    		[ new_id, it[1] ]
    	}
    	demulti_fastqs = bam2Fastq(demulti_bams)

	emit:
    	basecalled_fastq = demulti_fastqs.groupTuple()
    	demulti_report = dem_res.bar_summary
 
}

 workflow BASECALL_DEMULTIMOD {
    take: 
    input_fast5
    model_folder
    
    main:
     	model_folders = downloadModel(input_fast5.first().combine(model_folder))
        models = model_folders.collect().map{ [ it ] }
    	bam = baseCall(input_fast5.combine(models))
    	dem_res = demultiPlex(bam)
    	demulti_bams = dem_res.demulti_bams.transpose().map{
            def bam_name = it[1].getSimpleName()
            def barcode = "${bam_name}".split("_").last()
    	    def new_id = "${it[0]}.${barcode}"
    		[ new_id, it[1] ]
    	}
    	demulti_fastqs = bam2ModFastq(demulti_bams)

	emit:
    	basecalled_fastq = demulti_fastqs.groupTuple()
    	demulti_report = dem_res.bar_summary
 
}


 workflow BASECALL {
    take: 
    input_fast5
    model_folder
    
    main:
    	model_folders = downloadModel(input_fast5.first().combine(model_folder))
        models = model_folders.collect().map{ [ it ] }
    	bam = baseCall(input_fast5.combine(models)).basecalled_bam
    	bam2Fastq(bam)

	emit:
    	basecalled_fastq = bam2Fastq.out.basecalled_fastq
  
}




 workflow BASECALLMOD {
    take: 
    input_fast5
    model_folder
    
    main:
    	model_folders = downloadModel(input_fast5.first().combine(model_folder))
        models = model_folders.collect().map{ [ it ] }
    	bam = baseCall(input_fast5.combine(models)).basecalled_bam
    	bam2ModFastq(bam)

	emit:
    	basecalled_fastq = bam2ModFastq.out.basecalled_fastq
 
}





workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
} 


