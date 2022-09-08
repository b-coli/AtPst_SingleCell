#!/bin/bash

## Read in input arguments
sample=$1

## Run velocyto using default parameters fr 10x and edited gtf
singularity exec docker://dynverse/velocyto:latest velocyto run10x -@ 15 data/expression/${sample}/ genomes/compiled_references/at10/at10.gtf

