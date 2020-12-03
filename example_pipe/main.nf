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
outputMap = "${params.output}/Aln_file"

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
	include { GET_VERSION as FQ_VER; FASTQCP } from "${subwdir}/qc/fastqc" addParams(OUTPUT: outputQC) 
	include { GET_VERSION as BWA_VER; BWA_ALL } from "${subwdir}/alignment/bwa" addParams(OUTPUT: outputMap, LABEL: 'big_mem_cpus') 
	aln = BWA_ALL(params.reference, fastq_files)
	all_ver = FQ_VER().mix(BWA_VER()).collectFile(name: 'tool_version.txt', newLine: false, storeDir:outputQC)
}

