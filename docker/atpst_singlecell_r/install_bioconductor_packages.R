BiocManager::install(
  c(
    'batchelor', 
    'BiocGenerics', 
    'DelayedArray', 
    'DelayedMatrixStats',
    'DESeq2',
    'DropletUtils', 
    'edgeR', 
    'ggrastr',
    'HDF5Array', 
    'limma', 
    'lme4', 
    'Matrix.utils',
    'S4Vectors', 
    'SingleCellExperiment',
    'SummarizedExperiment', 
    'terra',
    'topGO'
  ),
  version = "3.12",
  update = FALSE)