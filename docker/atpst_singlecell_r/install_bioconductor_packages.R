BiocManager::install(
  c(
    'batchelor', 
    'BiocGenerics', 
    'DelayedArray', 
    'DelayedMatrixStats',                     
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
    'terra'
  ),
  version = "3.12",
  update = FALSE)