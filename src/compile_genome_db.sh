#!/usr/bin/sh

## Download Araport11 annotations and genome fasta from phytozome

## Set up folder
mkdir genomes/compiled_references/at10

## Unzip phytozome files
zcat genomes/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/assembly/Athaliana_447_TAIR10.fa.gz > genomes/compiled_references/at10/at10.fa
zcat genomes/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation/Athaliana_447_Araport11.gene_exons.gff3.gz > genomes/compiled_references/at10/at10.gff3

## Format annotation file, edit gtf files to contain gene id information needed
singularity exec docker://dceoy/gffread:latest gffread -g genomes/compiled_references/at10/at10.fa -T genomes/compiled_references/at10/at10.gff3 | \
	sed -E 's/gene_name [[:graph:]]+;//' | \
      	sed -E 's/.Araport11.447//g' > \
      	genomes/compiled_references/at10/at10.gtf
	
## Setup cellragner db
singularity exec docker://cumulusprod/cellranger:6.0.1 cellranger mkref \
  --genome=at10 \
  --fasta=genomes/compiled_references/at10/at10.fa \
  --genes=genomes/compiled_references/at10/at10.gtf

## Clean up
mv at10/* genomes/compiled_references/at10
mv Log.out log/cellranger_mkref_log.out
rmdir at10/
