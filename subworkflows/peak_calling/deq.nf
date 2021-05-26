/*
*  DEQ 
*/

 
params.LABEL = ""
params.EXTRAPARS = ""

params.OUTPUT = ""
params.CONTAINER = "biocorecrg/deq:0.1"

include { unzipCmd } from '../global_functions.nf'


process getVersion {
    container params.CONTAINER

    output:
	stdout emit: out    
    
    shell:
    """
    echo moaims' '`Rscript -e "library('deq'); packageVersion('deq')"` 2>/dev/null
    """
}


process de {
    label (params.LABEL)
    tag { comp_id }
    container params.CONTAINER
    if (params.OUTPUT != "") { publishDir(params.OUTPUT, mode:'copy') }

    input:
    tuple val(comp_id), path(peaks), path(input_untreated), path(ip_untreated), path(input_treated), path(ip_treated)
    path(annofile)
    val(read_size)
    val(frag_len)
    

    output:
    tuple val(comp_id), path("${comp_id}_deq_results.txt"), emit: deq_res
    
	script:
  def input_untreated_r = "c(\'${input_untreated.join('\',\'')}\')"
  def ip_untreated_r = "c(\'${ip_untreated.join('\',\'')}\')"
  def input_treated_r = "c(\'${input_treated.join('\',\'')}\')"
  def ip_treated_r = "c(\'${ip_treated.join('\',\'')}\')"
  def outfile = "\'${comp_id}_deq_results.txt\'"
  def gtffile = "\'${annofile}\'"
  def peaksfile = "\'${peaks}\'"
  def cmd = "input.bams, ip.bams, treated.input.bams, treated.ip.bams, peak.files, gtf, outfi = out, readlen = ${read_size}, fraglen = ${frag_len}, nthreads = ${task.cpus}"
  if (params.EXTRAPARS != "") {
    cmd = "${cmd}, ${params.EXTRAPARS}"
  }
     """ 
    Rscript -e "library('deq') 
    input.bams <- ${input_untreated_r}
    ip.bams <- ${ip_untreated_r}
    treated.input.bams <- ${input_treated_r}
    treated.ip.bams <- ${ip_treated_r}
    peak.files <- ${peaksfile}
    gtf <- ${gtffile}
    out <- ${outfile}
    deq(${cmd})"

     """
}



workflow DE {
    take: 
    data_de
    annotation
    read_size
    frag_len
    
    main:
		out = de(data_de, annotation, read_size, frag_len)
    emit:
    	out
}




workflow GET_VERSION {
    main:
		getVersion()
    emit:
    	getVersion.out
}    
