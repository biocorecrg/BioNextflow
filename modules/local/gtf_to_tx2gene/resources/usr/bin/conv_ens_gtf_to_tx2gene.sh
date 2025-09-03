#!/bin/bash

if [ x"$1" == x ]; then

        echo "please specify a GTF file"

        exit 1

fi

if [[ $1 =~ ".gz" ]]; then   
    catcmd="zcat"

else 
    catcmd="cat"

fi

$catcmd $1 |\
awk -F"\t" ' {\
	if ($3=="exon") {split($9,array,";"); \
		for (i in array) { \
			if (array[i]~"gene_id") dict["1"]=array[i]; \
			if (array[i]~"transcript_id") dict["2"]=array[i]; \
			gsub(/gene_id\ /, "", dict["1"]);\
			gsub(/transcript_id\ /, "", dict["2"]);\
		}\
		print dict[2]"\t"dict[1];
	} \
}' | sed s/\"//g| sed s/\ //g | sort | uniq > tx2gene.csv
