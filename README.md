# ![BioNextflow](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) BioNextflow

Groovy classes to be easily embedded in nextflow worflows. 
The idea is to make a repository of classes that groups methods for simplifying the portability of the code among different pipelines. 

In brief, insead of using the same "hard-coded" command line in different pipelines you can call the class.method().
You need to install the library in your nextflow pipeline for being automatically included in your nextflow script.


### Code before

    """
    if [ `echo ${reference_file} | grep ".gz"` ]; then 
       zcat ${reference_file} > ${index}.fa			
       bowtie-build --threads ${task.cpus} ${index}.fa ${index}
       else 
            ln -s ${reference_file} ${index}.fa 
            bowtie-build --threads ${tasks.cpus} ${index}.fa ${index}
       fi
    """

### Code after
    def aligner = new NGSaligner(reference_file:genome_file, index:"genome_index", cpus:task.cpus)
    aligner.doIndexing("bowtie")

### Installation
Edit the INSTALL.sh file to select the required version

    sh INSTALL.sh 

### Get started 

To run the test execute the following command: 

    nextflow run testFunctions.nf
