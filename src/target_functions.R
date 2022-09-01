read_bulk_dges <- function(sns) {
  dge_list <- map(sns, read_bulk_dge) %>%
    reduce(left_join, by = "Locus") %>% 
    column_to_rownames("Locus") %>%
    .[,sns]
}

read_bulk_dge <- function(sn) {
  fn <- paste0("data/expression/", sample, "/", sample, "ReadsPerGene.out.tab")
  
  expression <- read_tsv(
    filename, 
    skip = 4, 
    col_names = c("Locus", "expr", sample, "expr_stranded_rev"), 
    col_types = cols())
  
  expression %>% dplyr::select(Locus, one_of(sample))
}