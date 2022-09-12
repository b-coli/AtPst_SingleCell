library(targets)
library(tarchetypes)
library(edgeR)
library(Seurat)
library(tidyverse)

source("src/target_functions.R")

sample_metadata <- readr::read_csv("metadata/sample_metadata.csv")

bulk_data <- list(
  tar_target(
    bulk_dge, read_bulk_dges(
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
  
  tar_target(
    protoplast_loci,
    bulk_de_table %>%
      filter(de_type_proto == "Proto Down" | 
             de_type_proto == "Proto Up") %>% 
      pull(Locus)
  ),
  
  tar_render(
    bulk_analysis_report, 
    path = "src/bulk_analysis_report.Rmd", 
    output_dir = "reports/")
)

individual_objects <- list(
  tar_target(
    leaf_ref, 
    build_reference_dataset(path = "data/expression/KimPublished")
  ),
  
  tar_map(
    values = filter(sample_metadata, Sequence_Type == "scRNAseq"),
    names = "Sample_Name",
    unlist = F,
    tar_target(
      loom_file,
      paste0("data/expression/", Sample_Name, "/velocyto/", Sample_Name, ".loom"),
      format = "file"
    ),
    
    tar_target(
      raw_sobj,
      loom_to_sobj(
        loom_file = loom_file,
        md_file = Cell_Metadata_File
      )
    ),
    
    tar_target(sobj, process_sobj(sobj = raw_sobj, ref = leaf_ref))
  )
)

integrated_objects <- list(
  tar_target(
    atpst_sobj,
    integrate_sobjs(
      sobj_list = list(sobj_Control, sobj_DC3000),
      exclude_organelle_loci = T,
      exclude_features = protoplast_loci
    )
  ),
  
  tar_target(
    leaf_sobj, 
    integrate_sobjs(
      sobj_list = list(
        sobj_Control, 
        sobj_DC3000,
        sobj_KimLeaf1,
        sobj_KimLeaf2,
        sobj_ProckoLeaf1,
        sobj_ProckoLeaf2,
        sobj_Lopez.AnidoLeaf1),
      exclude_organelle_loci = T,
      exclude_features = protoplast_loci,
      method = "rpca"
    )
  )
)

list(
  bulk_data,
  individual_objects,
  integrated_objects
)
