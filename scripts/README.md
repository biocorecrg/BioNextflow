Here some useful script:

- get_singularity.py
	get singularity images from a json file output from nextflow inspect


For using it, you need to run the pipeline as usual but instead of 'nextflow run' you have to use 'nextflow inspect'. Then you can run on the json file the get_singularity.py script.

```
nextflow inspect -params-file params.yaml -profile slurm main.nf  > inspect.json

python get_singularity.py -j inspect.json -c ./singularity_folder 
```

