/*
*  nanoplot module
*  This workflow allows to make QC on input data
*  It needs input fastq
*/

params.LABEL = ""
params.CONTAINER = "biocorecrg/mopprepr:0.8"
params.OUTPUT = ""

process MOP_nanoPlot {
    tag { id }
    label (params.LABEL)
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }
   
    input:
    tuple val(id), path(bamfile)

    output:
    tuple val(id), path("*_plot/*"), emit: plots optional true
    tuple val(id), path("${id}_stats_mqc.png"), emit: pngs optional true 
    
    script:
    """
    NanoPlot --bam ${bamfile} -o ${id}_plot --maxlength 5000 -t ${task.cpus}
    mkdir tmp_dir
    cp ${id}_plot/PercentIdentityvsAverageBaseQuality_kde.png tmp_dir
    cp ${id}_plot/LengthvsQualityScatterPlot_dot.png tmp_dir
    cp ${id}_plot/HistogramReadlength.png tmp_dir 
    cp ${id}_plot/Weighted_HistogramReadlength.png tmp_dir
    gm montage tmp_dir/*.png -tile 2x2 -geometry 800x800 ${id}_stats_mqc.png
    rm -fr tmp_dir
    """
}



workflow MOP_QC {
    take: 
    fastq
    
    main:
    out = MOP_nanoPlot(fastq)
    
    emit:
   	plots = out.plots
    pngs = out.pngs

}

process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    NanoPlot --version
    """
}


workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}   
 
