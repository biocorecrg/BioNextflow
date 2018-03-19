#!/usr/bin/env nextflow 
/*
 * Copyright (c) 2018, Centre for Genomic Regulation (CRG).
 *
 *   This file is part of 'BioNextflow'.
 *
 *   BioNextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   BioNextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with BioNextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 *
 * @authors
 * Luca Cozzuto <lucacozzuto@gmail.com>
 *
 * 
 */

/* 
 * Check that helper function `mappingPairsWithSTAR` returns the exepcted result
 */

assert BiologicalFunctions.mappingPairsWithSTAR("idA","genome.fa","fastq1",4) == """
        	if [ `echo no == "debug"` ]; then print="echo "; else print=""; fi
			if [ `echo fastq1 | grep ".gz"` ]; then gzip="--readFilesCommand zcat"
			else gzip=""
			fi
            	\$print STAR --genomeDir genome.fa \
                     --readFilesIn fastq1 \
                     \$gzip \
                     --outSAMunmapped None \
                     --outSAMtype BAM SortedByCoordinate \
                     --runThreadN 4 \
                     --quantMode GeneCounts \
                     --outFileNamePrefix idA

                \$print mkdir STAR_idA
                \$print mv idAAligned* STAR_idA/.
                \$print mv idASJ* STAR_idA/.
                \$print mv idAReadsPerGene* STAR_idA/.
                \$print mv idALog* STAR_idA/.   
        """
        .stripIndent()
