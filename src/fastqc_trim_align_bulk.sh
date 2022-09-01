#!/bin/bash

## Pre-processes bulk RNA-seq data. Assumes data is a single set of fastq files (i.e. sequenced across a single lane),
## so that for single-end runs, just a single fastq file is present in the fastq folder, for paired-end, just two files
## are present. Takes as a sole input argument the sample name, and expects to find fastq file(s) in a folder at:
## data/fastq/SAMPLE_NAME.

sample_name=$1

## Run FastQC on the raw fastq files
singularity exec docker://staphb/fastqc:0.11.9 fastqc \
	-t 15 \
	-o log/ \
	-f fastq \
	data/fastq/${sample_name}/*.fastq.gz

## Quality, length, and adapter trim using TrimGalore
singularity exec docker://dceoy/trim_galore:latest trim_galore -q 30 \
	-j 15 \
	-o data/fastq/${sample_name}/ \
	--gzip --stringency 4 \
	data/fastq/${sample_name}/*.fastq.gz

## Move trimming report to the log folder
mv data/fastq/${sample_name}/*report.txt log/

## Run FastQC on trimmed fastq files
singularity exec docker://staphb/fastqc:0.11.9 fastqc -o log/ -f fastq data/fastq/${sample_name}/*trimmed*.fq.gz

## Create the expression directory
mkdir -p data/expression/${sample_name}/

## Use the STAR aligned (in the cellranger container) to align reads
singularity exec docker://cumulusprod/cellranger:6.0.1 /software/cellranger-6.0.1/lib/bin/STAR \
  --genomeDir genomes/compiled_references/at10/star \
  --runThreadN 15 \
  --readFilesCommand zcat \
  --runMode alignReads \
  --quantMode GeneCounts \
  --limitOutSJcollapsed 4000000 \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 16000000000 \
  --readFilesIn data/fastq/${sample_name}/*trimmed*.fq.gz \
  --outFileNamePrefix data/expression/${sample_name}/${sample_name}_
