library(targets)
library(tarchetypes)
library(edgeR)
library(Seurat)
library(tidyverse)

source("src/target_functions.R")

sample_metadata <- readr::read_csv("metadata/sample_metadata.csv")

bulk_data <- list(
  tar_target(bulk_dge, 
             read_bulk_dges(sns = sample_metadata %>%
                              filter(Sequence_Type == "RNAseq") %>%
                              pull(Sample_Name)
                            )
             ),
  tar_target(bulk_de_table, run_edger(dge = bulk_dge, 
                                  md = sample_metadata %>% 
                                    filter(Sequence_Type == "RNAseq")))
  
)

list(
  bulk_data
)
