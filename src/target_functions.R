read_bulk_dges <- function(sns) {
  dge_list <- map(sns, read_bulk_dge) %>%
    reduce(left_join, by = "Locus") %>% 
    column_to_rownames("Locus") %>%
    .[,sns]
}

read_bulk_dge <- function(sn) {
  fn <- paste0("data/expression/", sn, "/", sn, "_ReadsPerGene.out.tab")
  
  expression <- read_tsv(
    fn, 
    skip = 4, 
    col_names = c("Locus", "expr", sn, "expr_stranded_rev"), 
    col_types = cols())
  
  expression %>% dplyr::select(Locus, one_of(sn))
}

run_edger <- function(dge, md) {
#  md <- md %>% 
#    select(Sample_Name, Treatment, Sample_Type) %>% 
#    group_by(Treatment, Sample_Type) %>% 
#    mutate(Replicate = 1:n())
  
  erobj_norm <- generate_er_object(dge, md)
  
  proto_idx <- which(str_detect(colnames(erobj_norm$design), "Proto"))
  treat_idx <- which(str_detect(colnames(erobj_norm$design), "DC3000"))
  
  de_genes_proto <- edger_glmqfit_wrapper(
    x = erobj_norm, 
    design = erobj_norm$design, 
    coef_idx = proto_idx,
    suffix = "proto"
  )
  
  de_genes_treat <- edger_glmqfit_wrapper(
    x = erobj_norm, 
    design = erobj_norm$design, 
    coef_idx = treat_idx,
    suffix = "treat"
  )
  
  de_table <- full_join(de_genes_proto, de_genes_treat) %>%
    mutate(de_type_treat = 
             case_when(logFC_treat < -2 & p_adj_treat < 0.01 ~ "DC3000 Down",
                       logFC_treat > 2 & p_adj_treat < 0.01 ~ "DC3000 Up",
                       abs(logFC_treat) <= 2 | p_adj_treat >= 0.01 ~ "Not DE"),
           de_type_proto = 
             case_when(logFC_proto < -0.5 & p_adj_proto < 0.01 ~ "Proto Down",
                       logFC_proto > 0.5 & p_adj_proto < 0.01 ~ "Proto Up",
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