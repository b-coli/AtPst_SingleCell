# AtPst_SingleCell
Single-cell analysis of pseudomonas syringae DC3000 infection of Arabidopsis thaliana leaf tissue.

## Prerequisites
Requires a unix/linux environment with docker

This analysis makes use of several dockerhub images:
  * dceoy/gffread:latest          GFFRead utility for converting GFF3 files to GTF format
  * dceoy/trim_galore:latest      Trimming program for bulk RNA-seq reads
  * cumulusprod/cellranger:6.0.1  CellRanger tool (10x Genomics) for aligning reads and assigning cells, contains STAR
  * bcoli/AtPst_SingleCell:1.0    Custom R environment for analyzing single-cell data
  * mparikhbroad/velocyto:1.0.1   Velocyto software for determining spliced/unspliced counts

