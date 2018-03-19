#!/usr/bin/env nextflow

/*
 * Copyright (c) 2018, Centre for Genomic Regulation (CRG) and the authors.
 *
*
 */

/* 
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
