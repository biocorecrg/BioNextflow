# ![BioNextflow](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) BioNextflow

BioNextflow is a collection of sub-workflows that can be used in any [Nextflow DSL2 pipeline](https://www.nextflow.io/docs/latest/dsl2.html). They are created with the idea in mind of having a single file with different sub-workflow per tool, or combination of tools for a particular task. 

Ideally we will have a single file per tool, in which some custom parameters will indicate the containers, the command line extra parameters, the label etc. 

The fille will contain both processes and subworkflows that will be called from the main script. 

The input should be always a tuple val(id), path(inputfile) so that the id can be used for generating the output name if possible.

Luca

