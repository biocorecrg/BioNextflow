Here some useful script:

- get_singularity.py
	get singularity images from a json file output from nextflow inspect

```
nextflow inspect -params-file params.yaml -profile slurm main.nf  > inspect.json

python get_singularity.py -j inspect.json -c ./singularity_folder 
```

