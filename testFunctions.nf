#!/usr/bin/env nextflow

/*
 * Copyright (c) 2018, Centre for Genomic Regulation (CRG) and the authors.
 *
*
 */

/* 
 * Main ChIPSeq pipeline script for Bioinformatics Core @ CRG
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 *
 * 
 */

/*
 * 
 */

	process testFunction {
	echo true
	
	script:
		BiologicalFunctions.mappingPairsWithSTAR("idA","genome.fa","fastq1",4)
	
	}