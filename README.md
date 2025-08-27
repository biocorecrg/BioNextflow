<img align="right" href="https://biocore.crg.eu/" src="https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png" />

# ![BioNextflow](https://github.com/biocorecrg/BioNextflow/blob/master/bionextflow.png)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6394333.svg)](https://doi.org/10.5281/zenodo.6394333)
[![Nextflow version](https://img.shields.io/badge/Nextflow-21.04.1-brightgreen)](https://www.nextflow.io/)
[![Nextflow DSL2](https://img.shields.io/badge/Nextflow-DSL2-brightgreen)](https://www.nextflow.io/)


BioNextflow is a collection of modules and sub-workflows that can be used in any [Nextflow DSL2 pipeline](https://www.nextflow.io/docs/latest/dsl2.html). They are created with the idea in mind of having a single file per tool, containing different sub-workflows and modules.

The version v3.0 follows the rules from nf-core (less strict).

General functions are placed in **global_functions.nf** file and can be called from modules / subworkflows.

# Installation
The modules can be downloaded and used as they are or can be embedded as submodules.

## Submodules
To embed the submodules you can do:

```bash
git submodule add -b v3.0 https://github.com/biocorecrg/BioNextflow
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


