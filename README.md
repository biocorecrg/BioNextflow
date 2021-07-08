# ![BioNextflow](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) BioNextflow

BioNextflow is a collection of modules and sub-workflows that can be used in any [Nextflow DSL2 pipeline](https://www.nextflow.io/docs/latest/dsl2.html). They are created with the idea in mind of having a single file per tool, containing different sub-workflows and modules.

Ideally in each file we must define some custom parameters that indicate:
- the container 
- the labels
- extra parameters for the command line

The file will contain both processes and subworkflows that can be called from the main script. 

The input should be always a tuple with the id and a list such as:

```
tuple [val(id), [ val (something) ]
``` 

The id must be used for generating the output name if possible like:

```
${id}.out
```

Each final subworkflow should be able to handle both single end and paired end in a transparent way. This can be achieved by using as input this definition of single end:

```
tuple [val(id), [path(singleEND) ]
``` 

and this for paired ends:

```
tuple [val(id), [path(readA), path(readB) ]
```

Ideally there should be one or two separate modules for handling single reads and paired ends, while the final workflow should be smart enough to choose which module has to be used.

The input files should be always zipped. If unzipping is required the unzipped files should be removed after the end of the command execution.

```
zcat myfile.txt.gz > myfile.txt
do something myfile.txt > myfile.out
rm myfile.txt
```

Modules and subworkflows can be called from other subworkflows. General functions are placed in **global_functions.nf** file and can be called from modules / subworkflows.

# Installation
The modules can be downloaded and used as they are or can be embedded as submodules.

## Submodules
To embed the submodules you can do:

```bash
git submodule add https://github.com/biocorecrg/BioNextflow
...

Cloning into 'BioNextflow'...
remote: Enumerating objects: 636, done.
remote: Counting objects: 100% (306/306), done.
...

```

then you can add the .gitmodule file that is generated:

```bash
git add .gitmodules
git commit -m "adding modules"
git push
```

## Upgrading the submodules

You can upgrade the submodules using this command

```bash
git submodule update --remote --merge

git add BioNextflow
git commit -m "upgrading bioNextflow"
git push
```


