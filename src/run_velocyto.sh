#!/bin/bash

## Read in input arguments
sample=$1

## Run velocyto using default parameters fr 10x and edited gtf
singularity exec docker://mparikhbroad/velocyto:1.0.1 velocyto run10x data/expression/${sample}/ genomes/compiled_references/at10/at10.gtf

