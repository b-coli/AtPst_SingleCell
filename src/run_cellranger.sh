#!/bin/bash

## Read in input arguments
sample_name=$1
dir=data/fastq/${sample_name}/

## Run CellRanger using default parameters (and expected cells = 10000)
singularity exec docker://cumulusprod/cellranger:6.0.1 cellranger count \
  --nosecondary \
  --expect-cells 10000 \
  --include-introns \
  --localcores 15 \
  --id ${sample_name} \
  --transcriptome genomes/compiled_references/at10 \
  --fastqs $dir \
  --sample $sample_name

## Clean up
mv ${sample_name} data/expression/
