read_bulk_dges <- function(sns) {
  ## Read in multiple DGEs from STAR output
  ## sns = sample name
  ## Calls read_bulk_dge to format file name and read data
  dge_list <- map(sns, read_bulk_dge) %>%
    purrr::reduce(left_join, by = "Locus") %>% 
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
                       abs(logFC_proto) <= 0.5 | p_adj_proto >= 0.05 ~ "Not DE"))
  
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

loom_to_sobj <- function(loom_file, md_file = NA, sample_name = NA) {
  if(!file.exists(loom_file)) return(NA)
  
  dge <- SeuratWrappers::ReadVelocity(loom_file)
  sobj <- as.Seurat(dge)
  
  sobj[["RNA"]] <- sobj@assays$spliced
  DefaultAssay(sobj) <- "RNA"
  sobj$Sample_Name = sample_name
  
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

process_sobj <- function(sobj, ref, de_table) {
  if(is.na(sobj)) return(NA)
  cells <- sobj@meta.data %>%
    as_tibble(rownames = "Cell") %>%
    filter(nUMI >= UMI_Thresh) %>%
    filter(pct.mito < 0.01) %>%
    pull(Cell)
  
  sobj <- sobj[,cells]
  
  sobj <- seurat_pipeline(sobj, norm.method = "SCT")
  sobj <- transfer_labels(query = sobj, reference = ref, column = "Cell_Type")
  
  dc3000_up_genes <- filter(de_table, de_type_treat == "DC3000 Up") %>% pull(Locus)
  dc3000_down_genes <- filter(de_table, de_type_treat == "DC3000 Down") %>% pull(Locus)
  sobj <- Seurat::AddModuleScore(sobj, name = "DC3000.Up", features = list(dc3000_up_genes), assay = "RNA")
  sobj <- Seurat::AddModuleScore(sobj, name = "DC3000.Down", features = list(dc3000_down_genes), assay = "RNA")
  
  return(sobj)
  
}

seurat_pipeline <- function(sobj, norm.method = "SCT", num_pcs = 50, features = NULL) {
  if(norm.method == "SCT") {
    sobj <- Seurat::SCTransform(sobj)
  } else if(norm.method == "LN") {
    sobj <- Seruat::NormalizeData(sobj)
    sobj <- Seurat::ScaleData(sobj)
  }
  
  sobj <- RunPCA(sobj, npcs = num_pcs, features = features)
  sobj <- FindNeighbors(sobj, dims = 1:num_pcs)
  sobj <- FindClusters(sobj, resolution = 0.8)
  sobj <- FindClusters(sobj, resolution = 1.2)
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
  
  features <- map(sobj_list, row.names) %>% purrr::reduce(intersect)
  
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
    num_pcs = 50,
    features = features
  )
  
  sobj_integrate <- RunUMAP(
    sobj_integrated, 
    dims = 1:50,
    n.components = 3,
    reduction.key = "UMAP3D_", 
    reduction.name = "umap3d"
  )
  
  return(sobj_integrated)
}

monocle_pipeline <- function(sobj, cells, features, gene_ids) {
  
  gene_names <- gene_ids %>%
    group_by(Locus) %>%
    summarize(Gene = Gene[1])
  
  umap <- Seurat::Embeddings(sobj, reduction = "umap") %>% as_tibble(rownames = "Cell")
  
  cell_data <- sobj@meta.data %>%
    as_tibble(rownames = "Cell") %>%
    left_join(umap, by = "Cell") %>%
    mutate(seurat_clusters = as.numeric(as.character(seurat_clusters)))
  
  dge <- Seurat::GetAssayData(sobj, assay = "SCT", slot = "data")
  dge <- dge[features,]
  dge <- dge[,cells]
  
  md <- cell_data %>% 
    filter(Cell %in% cells) %>%
    column_to_rownames("Cell")
  
  cds <- monocle3::new_cell_data_set(
    dge,
    cell_metadata = md[cells,],
    gene_metadata = data.frame(
      Locus = features, 
      gene_short_name = features, 
      row.names = features, 
      stringsAsFactors = F)
    )
  
  cds <- monocle3::preprocess_cds(cds, num_dim = 5, norm_method = "none")
  cds <- monocle3::reduce_dimension(
    cds, reduction_method = "UMAP", 
    preprocess_method = "PCA", 
    umap.metric = "correlation", 
    umap.min_dist = 0.01
  )
  
  cds <- monocle3::cluster_cells(cds)
  
  start_cluster <- cell_data %>%
    mutate(diff_score = DC3000.Up1 - DC3000.Down1) %>%
    group_by(seurat_clusters) %>%
    summarize(mean_score = mean(diff_score)) %>%
    slice_min(mean_score, n = 1) %>%
    pull(seurat_clusters)
  
  start_cell <- cell_data %>%
    filter(seurat_clusters == start_cluster) %>% 
    filter(Sample_Name == "DC3000") %>%
    mutate(Signiture = (DC3000.Up1 - DC3000.Down1)) %>% 
    slice_min(Signiture, n = 1) %>%
    pull(Cell) %>% .[1]
  
  cds <- monocle3::learn_graph(
    cds, 
    learn_graph_control = list(
      minimal_branch_len = 15), 
    use_partition = F
  )
  
  cds <- monocle3::order_cells(cds, root_cells = start_cell)
  
  return(cds)
}

monocle_get_var_genes <- function(cds) {
  pseudotime_genes <- graph_test(
    cds, 
    neighbor_graph="principal_graph", 
    cores=4
  )
  
  return(pseudotime_genes)

}

download_and_read_gz <- function(url, type = "rds") {
  tmpfn <- paste0("tmp.gz")
  system("gunzip ", tmpfn)
  tmpfn <- "tmp"
  download.file(url, destfile = tmpfn)
  if(type == "csv") out <- read_csv(tmpfn)
  if(type %in% c("txt", "tsv")) out <- read_tsv(tmpfn)
  if(type == "rds" & compression == ".gz") out <- read_rds(tmpfn)

  file.remove(tmpfn)
  return(out)
}

extract_goi_pt_data <- function(sobj, pt, goi, nclusters = NULL) {
  goi_expr <- GetAssayData(sobj, assay = "SCT", slot = "data")[goi,] %>%
    as_tibble(rownames = "Locus") %>%
    pivot_longer(names_to = "Cell", values_to = "Expression", -Locus) %>%
    left_join(pt, by = "Cell") %>%
    filter(!is.na(Pseudotime)) %>%
    group_by(Locus) %>%
    filter(var(Expression) > 0) %>%
    mutate(Scaled_Expression = as.numeric(scale(Expression, center = T, scale = T))) %>%
    ungroup()
  
  goi_expr_mat <- goi_expr %>% 
    arrange(Pseudotime) %>%
    select(Cell, Locus, Scaled_Expression) %>%
    pivot_wider(names_from = Cell, values_from = Scaled_Expression) %>%
    column_to_rownames("Locus") %>%
    as.matrix()
  
  goi_expr_hclust <- dist(goi_expr_mat) %>% hclust()
  
  cell_order <- goi_expr %>% 
    select(Cell, Pseudotime) %>% 
    unique() %>% 
    arrange(Pseudotime) %>% 
    pull(Cell)
  
  gene_order <- goi_expr_hclust$labels[goi_expr_hclust$order]
  
  goi_expr <- goi_expr %>%
    mutate(Locus = factor(Locus, levels = gene_order)) %>%
    mutate(Cell = factor(Cell, levels = cell_order))
  
  if(!is.null(nclusters)) {
    clusters <- cutree(goi_expr_hclust, k = nclusters) %>% 
      enframe("Locus", "Module")
    goi_expr <- left_join(goi_expr, clusters)
  }
  
  return(goi_expr)
}

plot_pt_heatmap <- function(sobj, pt, goi, gene_names, nclusters = NULL) {
  goi_expr <- extract_goi_pt_data(sobj, pt, goi$Locus, nclusters = nclusters)
  goi <- mutate(goi, Locus = factor(Locus, levels = levels(goi_expr$Locus)))
  goi_expr <- left_join(goi_expr, goi)
  
  gene_labels <- left_join(goi_expr, gene_names) %>%
    select(Locus, Gene) %>%
    mutate(Gene = if_else(is.na(Gene), Locus, Gene)) %>%
    arrange(Locus) %>%
    deframe()
  
  plt <- goi_expr %>% 
    filter(!is.na(Pseudotime)) %>%
    mutate(Scaled_Expr = scales::squish(Scaled_Expression, range = c(-3,3))) %>%
    ggplot(aes(x = Cell, y = Locus)) + 
    geom_tile(aes(fill = Scaled_Expr)) +
    scale_fill_viridis_c(option = "plasma") +
    scale_y_discrete(labels = gene_labels) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  plt
}

get_cell_data <- function(sobj) {
  md <- sobj@meta.data %>% as_tibble(rownames = "Cell")
  reductions <- map(names(sobj@reductions), function(red) {
    Embeddings(sobj, reduction = red) %>%
      as_tibble(rownames = "Cell")}) %>%
    purrr::reduce(full_join)
  
  x = left_join(md, reductions)
  return(x)
}

read_archived_sobj <- function(fn, cluster_names) {
  sobj <- read_rds(fn)
  cell_data <- get_cell_data(sobj)
  cluster_renamed <- cluster_names %>%
    mutate(seurat_clusters = factor(
      seurat_clusters, 
      levels = levels(sobj@meta.data$seurat_clusters))) %>%
    right_join(cell_data) %>%
    select(Cell, Cluster) %>%
    column_to_rownames("Cell")
  
  sobj <- AddMetaData(sobj, cluster_renamed)
  
  return(sobj)
  
}