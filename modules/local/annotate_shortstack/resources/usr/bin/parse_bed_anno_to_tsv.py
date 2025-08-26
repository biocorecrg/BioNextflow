#!/usr/bin/env python

import csv
import gzip
import argparse
from collections import defaultdict

def parse_bed_to_tsv(bed_gz_file, tsv_file):
    """Parses a gzipped bed file and converts it into a tab-separated text file."""

    with gzip.open(bed_gz_file, 'rt') as bed, open(tsv_file, 'w', newline='') as tsv:
        writer = csv.writer(tsv, delimiter='\t', lineterminator='\n')
        annotation_headers = ["Chr", "Start", "End", "Name", "Score", "Strand",
                              "gene.id", "gene.name", "gene.type"]
        annotations   = {}
        writer.writerow(annotation_headers)  # Keep sample columns
        gene_names        = defaultdict(set)
        gene_bios         = defaultdict(set)
        exon_matches      = defaultdict(set)
        genes_matches     = defaultdict(set)

        for line in bed:

            # Process BED entry
            fields = line.strip().split("\t")
            chrom, start, end, name, score, strand = fields[:6]
            annotations[name] = "\t".join(fields[:6])
            # Kind of GTF anno that is overlapping
            annotype = fields[8]
            # Extract info field
            info = fields[-1]

            # default gene_name that can be missing
            gene_name     = "empty"
            gene_id       = "empty"
            gene_bio      = "empty"

            for entry in info.split(";"):
                entry = entry.lstrip()
                if entry.startswith("gene_id"):
                    gene_id = entry.split(" ")[1].replace('"', '')
                if entry.startswith("biotype") or entry.startswith("gene_type") or entry.startswith("gene_biotype"):
                    gene_bio = entry.split(" ")[1].replace('"', '')
                if entry.startswith("gene_name") or entry.startswith("Name"):
                    gene_name = entry.split(" ")[1].replace('"', '')

            if (annotype == "exon"):
                exon_matches[name].add(gene_id)

            if (annotype == "gene"):
                if (gene_name == "empty"):
                    gene_name = gene_id

                genes_matches[name].add(gene_id)
                gene_names[gene_id] = gene_name
                gene_bios[gene_id] = gene_bio


        for ann_name in annotations:

            my_exons    = exon_matches[ann_name]
            my_genes    = genes_matches[ann_name]
            my_introns  =  set(my_genes - my_exons)

            my_exnames  = {k: gene_names[k] for k in my_exons if k in gene_names}.values()
            my_intnames =  {k: gene_names[k] for k in my_introns if k in gene_names}.values()
            my_exbio    = {k: gene_bios[k] for k in my_exons if k in gene_bios}.values()
            my_intbio   = {k: "intron-" + gene_bios[k] for k in my_introns if k in gene_bios}.values()

            if (len(my_exons) > 0 and len(my_introns) > 0):
                id_list   = ",".join(my_exons) + "," + ",".join(my_introns)
                name_list = ",".join(my_exnames) + "," + ",".join(my_intnames)
                bio_list  = ",".join(my_exbio) + "," + ",".join(my_intbio)

            elif (len(my_exons) > 0):
                id_list   = ",".join(my_exons)
                name_list = ",".join(my_exnames)
                bio_list  = ",".join(my_exbio)

            elif (len(my_introns) > 0):
                id_list   = ",".join(my_introns)
                name_list = ",".join(my_intnames)
                bio_list  = ",".join(my_intbio)

            else:
                id_list   = "intergenic"
                name_list = "intergenic"
                bio_list  = "intergenic"

            new_fields = (annotations[ann_name] + "\t" + id_list + "\t" + name_list + "\t" + bio_list).split("\t")
            writer.writerow(new_fields)  # Keep sample columns

def main():
    parser = argparse.ArgumentParser(description="Convert a BED.gz file into a tab-separated text file.")
    parser.add_argument("-i", "--input", required=True, help="Input BED.gz file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")

    args = parser.parse_args()
    parse_bed_to_tsv(args.input, args.output)

if __name__ == "__main__":
    main()
