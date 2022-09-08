library(targets)
library(tarchetypes)
library(edgeR)
library(Seurat)
library(tidyverse)

source("src/target_functions.R")

sample_metadata <- readr::read_csv("metadata/sample_metadata.csv")

bulk_data <- list(
  tar_target(
    bulk_dge, 
    read_bulk_dges(
      sns = sample_metadata %>%
        filter(Sequence_Type == "RNAseq") %>%
        pull(Sample_Name)
    )
  ),
  tar_target(
    bulk_erobj_norm, 
    generate_er_object(
      dge = bulk_dge, 
      md = filter(sample_metadata, Sequence_Type == "RNAseq")
    )
  ),
  
  tar_target(
    bulk_de_table,
    get_de_genes(bulk_erobj_norm)
  ),
  
  tar_render(
    bulk_analysis_report, 
    path = "src/bulk_analysis_report.Rmd", 
    output_dir = "reports/", 
    knit_root_dir = getwd(),
    output_format = "github_document")
  
)


list(
  bulk_data
)
