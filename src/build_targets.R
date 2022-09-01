library(targets)
library(tarchetypes)
library(DESeq2)
library(Seurat)
library(tidyverse)

source("src/target_functions.R")

sample_metadata <- read_csv("metadata/sample_metadata.csv")

bulk_data <- list(
  tar_target(bulk_dge, 
             read_bulk_dges(sns = sample_metadata %>%
                              filter(Sequence_Type == "RNAseq") %>%
                              pull(Sample_Name)
                            )
             ),
  
)

list(
  bulk_data
)