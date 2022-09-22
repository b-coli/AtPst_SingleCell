library(targets)
library(tarchetypes)
library(edgeR)
library(Seurat)
library(tidyverse)
library(future)
library(monocle3)

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
      filter(de_type_proto %in% c("Proto Up", "Proto Down")) %>% 
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
        md_file = Cell_Metadata_File,
        sample_name = Sample_Name
      )
    ),
    
    tar_target(
      sobj, 
      process_sobj(
        sobj = raw_sobj, 
        ref = leaf_ref,
        de_table = bulk_de_table
      )
    )
  )
)

integrated_objects <- list(
  tar_target(
    atpst_sobj,
    integrate_sobjs(
      sobj_list = list(sobj_Control, 
                       sobj_DC3000),
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
      method = "cca"
    )
  )
)

pseudotime_objects <- list(
  tar_target(
    monocle_cells,
    archived_sobj@meta.data %>%
      as_tibble(rownames = "Cell") %>%
      filter(predicted.id == "Mesophyll") %>%
      filter(Sample_Name == "DC3000") %>%
      filter(seurat_clusters != 16) %>%
      pull(Cell)
  ),
  
  tar_target(
    monocle_genes,
    bulk_de_table %>%
      filter(de_type_proto != "Not DE") %>%
      filter(de_type_treat == "Not DE") %>%
      pull(Locus) %>%
      setdiff(x = row.names(archived_sobj@assays$SCT))
  ),
  
  tar_target(
    atpst_mobj,
    monocle_pipeline(
      sobj = archived_sobj,
      cells = monocle_cells,
      features = monocle_genes,
      gene_ids = gene_ids
    )
  ),
  
  tar_target(
    pt_var_genes,
    monocle_get_var_genes(atpst_mobj)
  )
)

archival_data <- list(
  tar_target(
    gene_ids_file,
    "data/gene_aliases_20200930.txt.gz",
    format = "file"
  ),

  tar_target(
    archived_sobj_file,
    "data/GSE213622_coaker_atpst_singlecell_sobj.rds",
    format = "file"
  ),

  tar_target(
    gene_ids,
    read_tsv(
      gene_ids_file,
      col_types = cols(),
      col_names = c("Locus", "Gene", "Full")
    )
  ),

  tar_target(
     archived_sobj,
     read_rds(archived_sobj_file)
  )

)

list(
  bulk_data,
  individual_objects,
  integrated_objects,
  archival_data,
  pseudotime_objects
)