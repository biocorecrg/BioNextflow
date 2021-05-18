# ![BioNextflow](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) BioNextflow

BioNextflow is a collection of sub-workflows that can be used in any [Nextflow DSL2 pipeline](https://www.nextflow.io/docs/latest/dsl2.html). They are created with the idea in mind of having a single file with different sub-workflow per tool, or combination of tools for a particular task. 

Ideally we will have a single file per tool, in which some custom parameters will indicate the containers, the command line extra parameters, the label etc. 

The file will contain both processes and subworkflows that will be called from the main script. 

The input should be always:

```tuple [val(id), [path(inputfile) ]``` 

or

```tuple [val(id), [path(readA), path(readB) ]```

and ideally there should be one or two separate modules for handling single reads and paired ends, while the workflow should be smart enough to choose the modules for SE or PE.

The id must be used for generating the output name if possible.

The input files should be always zipped. If unzipping is required the unzipped files should be removed after the end of the command execution


