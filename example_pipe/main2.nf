#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/* 
 * Define the pipeline parameters
 */


params.help            = false

log.info """
====================================================
Example pipeline for using sub-workflows
====================================================
reference                 : ${params.reference}
fastq                     : ${params.fastq}
output                    : ${params.output}
"""

// Help 
if (params.help) exit 1

// Output folders
outputQC  = "${params.output}/QC_files"
outputBWA = "${params.output}/BWA_file"
outputSAM = "${params.output}/SAMB_file"

subwdir   = "${baseDir}/../subworkflows"
/*
 * Creates the channels that emits fastq files
 * it can be paired ends or single ends
 */
Channel
 .fromFilePairs( params.fastq , size: ("${params.fastq}" =~ /\{/) ? 2 : 1)
 .ifEmpty { error "Cannot find any reads matching: ${params.fastq}" }
 .set { fastq_files }

workflow {	
	include { ALL as SALMON_ALL } from "${subwdir}/alignment/salmon" addParams(OUTPUT: outputBWA)
        //include { GET_VERSION as FQ_VER; FASTQCP } from "${subwdir}/qc/fastqc" addParams(OUTPUT: outputQC) 
	//include { GET_VERSION as BWA_VER; BWA_ALL } from "${subwdir}/alignment/bwa" addParams(OUTPUT: outputBWA, LABEL: 'big_mem_cpus') 
	//include { GET_VERSION as SAM_VER; SAMBLASTER_ALL } from "${subwdir}/alignment/samblaster" addParams(OUTPUT: outputSAM, LABEL: 'big_mem_cpus') 
	//include { FILTER } from "${subwdir}/trimming/trimmomatic" addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36") 

        SALMON_ALL(params.reference, fastq_files)
        //FILTER(fastq_files)
	//FASTQCP(fastq_files)
	//BWA_ALL(params.reference, fastq_files)
        //SAMBLASTER_ALL(params.reference, fastq_files)
	//FQ_VER().mix(BWA_VER(), SAM_VER()).collectFile(name: 'tool_version.txt', newLine: false, storeDir:outputQC)
}

workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
