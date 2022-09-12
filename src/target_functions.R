read_bulk_dges <- function(sns) {
  ## Read in multiple DGEs from STAR output
  ## sns = sample name
  ## Calls read_bulk_dge to format file name and read data
  dge_list <- map(sns, read_bulk_dge) %>%
    reduce(left_join, by = "Locus") %>% 
    column_to_rownames("Locus") %>%
    .[,sns]
}

read_bulk_dge <- function(sn) {
  ## Read a single digital gene expression (DGE) matrix from STAR output
  
  ## Define file name from the sample name
  fn <- paste0("data/expression/", sn, "/", sn, "_ReadsPerGene.out.tab")
  
  ## Read in the counts data
  expression <- read_tsv(
    fn, 
    skip = 4, 
    col_names = c("Locus", "expr", sn, "expr_stranded_rev"), 
    col_types = cols())
  
  expression %>% dplyr::select(Locus, one_of(sn))
}

get_de_genes <- function(erobj_norm) {
  ## Run differential expression analysis for bulk and protoplasted control
  ## and DC3000 treatments.
  
  ## erobj_norm is the normalized edgeR object
  
  ## Define the column number in the design matrix where protoplasted or 
  ## treatment samples can be found
  proto_idx <- which(str_detect(colnames(erobj_norm$design), "Proto"))
  treat_idx <- which(str_detect(colnames(erobj_norm$design), "DC3000"))
  
  ## Generate differential expression table for protoplasting
  de_genes_proto <- edger_glmqfit_wrapper(
    x = erobj_norm, 
    design = erobj_norm$design, 
    coef_idx = proto_idx,
    suffix = "proto"
  )
  
  ## Generate differential expression table for DC3000 treatment
  de_genes_treat <- edger_glmqfit_wrapper(
    x = erobj_norm, 
    design = erobj_norm$design, 
    coef_idx = treat_idx,
    suffix = "treat"
  )
  
  ## Merge differential expression tables, define thresholds.
  de_table <- full_join(de_genes_proto, de_genes_treat) %>%
    mutate(de_type_treat = 
             case_when(logFC_treat < -2 & p_adj_treat < 0.01 ~ "DC3000 Down",
                       logFC_treat > 2 & p_adj_treat < 0.01 ~ "DC3000 Up",
                       abs(logFC_treat) <= 2 | p_adj_treat >= 0.01 ~ "Not DE"),
           de_type_proto = 
             case_when(logFC_proto < -0.5 & p_adj_proto < 0.05 ~ "Proto Down",
                       logFC_proto > 0.5 & p_adj_proto < 0.05 ~ "Proto Up",
                       abs(logFC_proto) <= 0.5 | p_adj_proto >= 0.01 ~ "Not DE"))
  
  return(de_table)
  
}

edger_glmqfit_wrapper <- function(x, design, coef_idx, suffix = NULL) {
  fit <- edgeR::glmQLFit(x, design) %>% 
    edgeR::glmQLFTest(coef = coef_idx)
  
  fit_de <- fit$table %>% 
    as_tibble(rownames = "Locus") %>%
    mutate(p_adj = p.adjust(PValue, method = "BH"))
  
  if(!is.null(suffix)) fit_de <- rename_with(fit_de, paste0, -Locus, "_", suffix)

  return(fit_de)
}

generate_er_object <- function(dge, md) {
  treatments <- md$Treatment
  sample_types <- md$Sample_Type
  design <- model.matrix(~ sample_types + treatments)
  
  dge <- dge[,md$Sample_Name]
  
  erobj_raw <- edgeR::DGEList(dge)
  genes_to_keep <- edgeR::filterByExpr(erobj_raw, design = design)
  
  erobj_norm <- erobj_raw[genes_to_keep, ] %>% 
    edgeR::calcNormFactors() %>%
    edgeR::estimateDisp(design)
}

build_reference_dataset <- function(path) {
  dge <- Seurat::Read10X(path)
  md <- get_mito_md(dge)
  
  sobj <- Seurat::CreateSeuratObject(dge, meta.data = md) %>%
    Seurat::SCTransform(vars.to.regress = "pct.mito") %>%
    Seurat::RunPCA(npcs = 50, assay = "SCT") %>%
    Seurat::FindNeighbors(dims = 1:50, nn.method = "rann") %>%
    Seurat::FindClusters(resolution = 0.8) %>%
    Seurat::RunUMAP(dims = 1:50, min.dist = 0.1, n.neighbors = 30)
  
  loci <- c("AT3G48740", "AT5G23660", "AT1G22710", "AT5G06850", "AT5G41920",
            "AT1G77990", "AT3G24140", "AT2G46070", "AT4G21750", "AT1G51500",
            "AT1G68530", "AT1G29910", "AT2G05100", "AT1G28230", "AT3G54420", 
            "AT2G45190")
  
  sct_data <- Seurat::GetAssayData(sobj, assay = "SCT", slot = "data")[loci,] %>%
    tibble::as_tibble(rownames = "Locus") %>%
    tidyr::gather("Cell", "Expression", -Locus) %>%
    left_join(as_tibble(sobj@meta.data, rownames = "Cell"), by = "Cell")
  
  cluster_annotation <- sct_data %>% 
    group_by(Locus, seurat_clusters) %>% 
    summarize(pct.expressed = mean(Expression > 0)) %>% 
    pivot_wider(names_from = Locus, values_from = pct.expressed) %>% 
    mutate(Cell_Type = case_when(
      AT3G48740 > 0.1 & AT5G23660 > 0.1 ~ "Phloem", # SWEET11 and 12
      AT1G77990 > 0.5 & AT5G41920 > 0.1 ~ "Bundle Sheath", # SULTR2;2, SCL23
      AT1G22710 > 0.2 & AT5G06850 > 0.1 ~ "Companion", # SUC2 & FTIP1
      AT3G24140 > 0.2 & AT2G46070 > 0.2 ~ "Guard", # FAMA & MAPK12
      AT4G21750 > 0.2 & AT1G51500 > 0.2 ~ "Epidermis", # ATML1 & CER5
      AT1G28230 > 0.1 & AT3G54420 > 0.1 ~ "Hydathode", # PUP1 & EP3
      AT2G05100 > 0.98 & AT1G29910 > 0.98 ~ "Mesophyll" # LCHB1.2, CAB3
    )) %>%
    select(seurat_clusters, Cell_Type) %>%
    left_join(as_tibble(sobj@meta.data, rownames = "Cell")) %>%
    select(Cell, Cell_Type) %>%
    column_to_rownames("Cell")
  
  sobj <- AddMetaData(sobj, cluster_annotation)
  sobj <- sobj[,!is.na(sobj$Cell_Type)]
  
  return(sobj)
}

get_mito_md <- function(dge) {
  mito_idx <- str_detect(row.names(dge), "ATMG.....")
  
  mito_sums <- Matrix::colSums(dge[mito_idx,]) %>%
    enframe("Cell", "nMt")
  
  all_sums <- Matrix::colSums(dge) %>%
    enframe("Cell", "nUMI")
  
  cell_md <- full_join(all_sums, mito_sums, by = "Cell") %>%
    mutate(pct.mito = nMt/nUMI) %>%
    column_to_rownames("Cell")
}

loom_to_sobj <- function(loom_file, md_file = NA) {
  if(!file.exists(loom_file)) return(NA)
  dges <- SeuratWrappers::ReadVelocity(loom_file)
  sobj <- CreateSeuratObject(dges$spliced)
  if(!is.na(md_file)) {
    md <- read_csv(md_file) %>% column_to_rownames("Barcode")
    sobj <- AddMetaData(sobj, md)
  }
  
  cell_stats <- get_cell_stats(sobj)
  sobj <- AddMetaData(sobj, cell_stats)
  
  return(sobj)
}

get_cell_stats <- function(sobj) {
  dge <- sobj@assays$RNA@counts
  all_loci <- row.names(dge)
  mito_loci <- all_loci[stringr::str_detect(all_loci, "ATMG.....")]
  chloro_loci <- all_loci[stringr::str_detect(all_loci, "ATCG.....")]
  
  umi_sums <- Matrix::colSums(dge) %>% tibble::enframe("Cell", "nUMI")
  mito_sums <- Matrix::colSums(dge[mito_loci,]) %>% tibble::enframe("Cell", "nMt")
  chloro_sums <- Matrix::colSums(dge[chloro_loci,]) %>% tibble::enframe("Cell", "nCp")
  
  cell_stats <- umi_sums %>%
    left_join(chloro_sums, by = "Cell") %>% 
    left_join(mito_sums, by = "Cell") %>%
    mutate(pct.mito = nMt/nUMI)
  
  umi_thresh <- cell_stats %>%
    filter(nUMI < 50000) %>%
    slice_max(nUMI, n = 10000*0.01) %>%
    pull(nUMI) %>%
    min()*0.1
  
  cell_stats$UMI_Thresh <- umi_thresh
  
  cell_stats <- column_to_rownames(cell_stats, "Cell")
  
  return(cell_stats)
}

process_sobj <- function(sobj, ref) {
  if(is.na(sobj)) return(NA)
  sobj <- sobj[,sobj$nUMI >= sobj$UMI_Thresh]
  sobj <- sobj[,sobj$nUMI < 50000]
  sobj <- sobj[,sobj$pct.mito < 0.01]
  
  sobj <- seurat_pipeline(sobj, norm.method = "SCT")
  sobj <- transfer_labels(query = sobj, reference = ref, column = "Cell_Type")
  
  return(sobj)
  
}

seurat_pipeline <- function(sobj, norm.method = "SCT", num_pcs = 50) {
  if(norm.method == "SCT") {
    sobj <- Seurat::SCTransform(sobj)
  } else if(norm.method == "LN") {
    sobj <- Seruat::NormalizeData(sobj)
    sobj <- Seurat::ScaleData(sobj)
  }
  
  sobj <- RunPCA(sobj, npcs = num_pcs)
  sobj <- FindNeighbors(sobj, dims = 1:num_pcs)
  sobj <- FindClusters(sobj, resolution = 0.8)
  sobj <- RunUMAP(sobj, dims = 1:num_pcs)
  
  return(sobj)
}

transfer_labels <- function(query, reference, column = "Cell_Type") {
  
  anchors <- FindTransferAnchors(reference = reference,
                                 query = query,
                                 normalization.method = "SCT")
  
  predictions <- TransferData(anchorset = anchors,
                              refdata = FetchData(reference, column) %>% as_tibble(rownames = "Cell") %>% deframe())
  
  sobj <- AddMetaData(query, metadata = predictions)
  
  return(sobj)
}

integrate_sobjs <- function(sobj_list, exclude_organelle_loci = TRUE, 
                            exclude_features = NULL, method = "cca") {
  
  features <- map(sobj_list, row.names) %>% reduce(intersect)
  
  if(!is.null(exclude_features)) features <- setdiff(features, exclude_features)
  
  if(exclude_organelle_loci) {
    organelle_loci <- str_detect(features, "AT[MC]G.....")
    features <- setdiff(features, organelle_loci)
  }

  sobj_list <- PrepSCTIntegration(sobj_list, anchor.features = features)
  
  sobj_integrated <- FindIntegrationAnchors(
    object.list = sobj_list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = method
  )
  
  sobj_integrated <- IntegrateData(
    anchorset = sobj_integrated, 
    normalization.method = "SCT"
  )
  
  sobj_integrated <- seurat_pipeline(
    sobj_integrated, 
    norm.method = "none", 
    num_pcs = 50
  )

  return(sobj_integrated)
}