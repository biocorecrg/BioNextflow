/*
* trimmomatic module
*/

params.CONTAINER = "quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4"
params.OUTPUT = "trimmomatic_output"

process trimmomatic {
    publishDir(params.OUTPUT, mode: 'copy')
    tag {pair_id }
    container params.CONTAINER

    input:
    tuple val (pair_id), path(fq1), path(fq2) 

    output:
    tuple val (pair_id), path("${pair_id}_1P.fq.gz"), path("${pair_id}_2P.fq.gz")

    script:

    """
    trimmomatic \
    PE -phred33 \
    $fq1 \
    $fq2 \
    ${pair_id}_1P.fq.gz \
    ${pair_id}_2P.fq.gz \
    ${pair_id}_1UP.fq.gz \
    ${pair_id}_2UP.fq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}
