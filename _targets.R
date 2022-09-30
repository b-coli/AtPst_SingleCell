library(targets)
library(tarchetypes)
library(edgeR)
library(Seurat)
library(future)
library(monocle3)
library(topGO)
library(tidyverse)

source("src/target_functions.R")
future::plan(strategy = "multisession", workers = 4)
options(future.globals.maxSize = 10*(1024^3))

sample_metadata <- readr::read_csv("metadata/sample_metadata.csv")

bulk_data <- list(
  tar_target(bulk_dge, read_bulk_dges(sns = sample_metadata %>% filter(Sequence_Type == "RNAseq") %>% pull(Sample_Name))),
  tar_target(bulk_erobj_norm, generate_er_object(dge = bulk_dge, md = filter(sample_metadata, Sequence_Type == "RNAseq"))),
  tar_target(bulk_de_table, get_de_genes(bulk_erobj_norm)),
  tar_target(protoplast_loci, bulk_de_table %>% filter(de_type_proto %in% c("Proto Up")) %>% pull(Locus)),
  tar_render(bulk_analysis_report, path = "src/bulk_analysis_report.Rmd", output_dir = "reports/")
)

individual_objects <- list(
  tar_target(leaf_ref, build_reference_dataset(path = "data/expression/KimPublished")),
  tar_map(
    values = filter(sample_metadata, Sequence_Type == "scRNAseq"),
    names = "Sample_Name",
    unlist = F,
    tar_target(loom_file, paste0("data/expression/", Sample_Name, "/velocyto/", Sample_Name, ".loom"), format = "file"),
    tar_target(raw_sobj, loom_to_sobj(loom_file = loom_file, md_file = Cell_Metadata_File, sample_name = Sample_Name)),
    tar_target(sobj, process_sobj(sobj = raw_sobj, ref = leaf_ref, de_table = bulk_de_table))
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
        sobj_Lopez.AnidoLeaf1,
        sobj_ZhangLeaf1),
      exclude_organelle_loci = T,
      exclude_features = protoplast_loci,
      method = "cca"
    )
  ),
  
  tar_target(cluster_deg, FindAllMarkers(archived_sobj, assay = "SCT", logfc.threshold = 0, min.pct = 0)),
  tar_target(cluster_sample_deg, find_markers2d(archived_sobj, comp_1 = "DC3000", comp_2 = "Control")),
  tar_render(pseudobulk_report, path = "src/pseudobulk_report.Rmd", output_dir = "reports/"),
  tar_render(dataset_stats_report, path = "src/dataset_stats_report.Rmd", output_dir = "reports/"),
  tar_render(combinedleaf_report, path = "src/combinedleaf_report.Rmd", output_dir = "reports/"),
  tar_render(de_report, path = "src/de_report.Rmd", output_dir = "reports/")
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
      filter(abs(logFC_proto) > 0.5 & p_adj_proto < 0.1) %>%
      filter(de_type_treat == "Not DE") %>%
      pull(Locus) %>%
      setdiff(x = row.names(archived_sobj@assays$SCT))
  ),
  
  tar_target(atpst_mobj, monocle_pipeline(sobj = archived_sobj, cells = monocle_cells, features = monocle_genes, gene_ids = gene_ids)),
  tar_target(pt_var_table, monocle_get_var_genes(atpst_mobj)),
  tar_target(pt_var_genes, pt_var_table %>% filter(morans_I > 0.2) %>% row.names()),
  tar_render(pseudotime_report, path = "src/pseudotime_report.Rmd", output_dir = "reports/")
)

archival_data <- list(
  tar_target(gene_ids_file, "data/gene_aliases_20200930.txt.gz", format = "file"),
  tar_target(gene_ids, read_tsv(gene_ids_file, col_types = cols(), col_names = c("Locus", "Gene", "Full"))),
  tar_target(archived_sobj_file, "data/GSE213622_coaker_atpst_singlecell_sobj.rds", format = "file"),
  tar_target(archived_sobj, read_archived_sobj(fn = archived_sobj_file, cluster_names = archived_clusterNames)),
  tar_target(archived_clusterNames_file, "data/cluster_renaming.csv", format = "file"),
  tar_target(archived_clusterNames, read_csv(archived_clusterNames_file))
)

go_information <- list(
  tar_target(
    gene_annotation, 
    download_and_read_gz(
      "https://arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt.gz", 
      type = "tsv",
      skip = 4,
      col_names = c(
        "Locus", "TAIR_Accession", "ObjectName", "Relationship", 
        "GO_Description", "GO", "TAIR_Keyword", "Aspect", "GOslim_Term", 
        "EvidenceCode", "EvidenceDescription", "SupportingEvidence", 
        "Reference", "Annotator", "Date"
      )
    )),
  tar_target(gene2GO, map_GO_to_genes(gene_annotation)),
  tar_target(marker_go, get_go_for_clusters(
    sobj = archived_sobj, 
    markers = cluster_deg, 
    min.pct = 0.1, 
    min.lfc = 1, 
    p = 0.01, 
    gene_to_go = gene2GO,
    exclude_loci = protoplast_loci)),
  tar_render(go_report, path = "src/go_report.Rmd", output_dir = "reports/")
  
)

list(
  bulk_data,
  individual_objects,
  integrated_objects,
  archival_data,
  pseudotime_objects,
  go_information
)
