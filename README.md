# ![BioNextflow](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) BioNextflow

Groovy classes to be easily embedded in nextflow worflows. 
The idea is to make a repository of classes that groups methods for simplifying the portability of the code among different pipelines. 

In brief, insead of using the same "hard-coded" command line in different pipelines you can call the class.method().
You need to install the library in your nextflow pipeline for being automatically included in your nextflow script.


### Code before

        """
        if [ `echo ${genome_file} | grep ".gz"` ]; then 
                zcat ${genome_file} > `basename ${genome_file} .gz`
                bowtie2-build --threads ${cpus} `basename ${genome_file} .gz` ${indexname}
                rm `basename ${genome_file} .gz`
                else bowtie2-build --threads ${cpus} ${genome_file} ${indexname}
                fi
        """

### Code after
    aligner = new(genome_file:genome_file, index:bowtie2genome", cpus:task.cpus)
    aligner.doAlignment()

### Installation
Edit the INSTALL.sh file to select the required version

  sh INSTALL.sh 

### Get started 

To run the test execute the following command: 

    nextflow run testFunctions.nf
