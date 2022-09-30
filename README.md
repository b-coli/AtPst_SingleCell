# AtPst_SingleCell
Single-cell analysis of pseudomonas syringae DC3000 infection of Arabidopsis thaliana leaf tissue.

## Prerequisites
Requires a unix/linux environment with docker

This analysis makes use of several dockerhub images:
|Tool                          |Purpose                                                                             | 
|------------------------------|------------------------------------------------------------------------------------|
|dceoy/gffread:latest          |GFFRead utility for converting GFF3 files to GTF format                             |
|dceoy/trim_galore:latest      |Trimming program for bulk RNA-seq reads                                             |
|cumulusprod/cellranger:6.0.1  |CellRanger tool (10x Genomics) for aligning reads and assigning cells, contains STAR|
|bcoli/AtPst_SingleCell:1.0    |Custom R environment for analyzing single-cell data                                 |                    
|dynverse/velocyto:latest      |Velocyto software for determining spliced/unspliced counts                          |
|ncbi/sra-tools                |NCBI SRA tools for downloading SRA data                                             |
|staphb/fastqc:0.11.9          |FastQC for evaluating fastq data quality                                            |

The analysis also uses resources on Phytozome, which will need to be downloaded from:

https://data.jgi.doe.gov/refine-download/phytozome?organism=Athaliana&expanded=447
 * Genome annotation: Athaliana_447_Araport11.gene_exons.gff3.gz
 * Genome sequence: Athaliana_447_TAIR10.fa.gz
 
And information from TAIR, which can be downloaded from:
https://www.arabidopsis.org/download_files/Public_Data_Releases/TAIR_Data_20200930/gene_aliases_20200930.txt.gz
https://arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt.gz
