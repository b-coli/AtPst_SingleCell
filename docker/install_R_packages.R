install.packages("https://cran.r-project.org/src/contrib/remotes_2.4.2.tar.gz")

remotes::install_version("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install(c(
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
  'pcaMethods', 
  'S4Vectors', 
  'SingleCellExperiment',
  'SummarizedExperiment', 
  'terra'),
  version = "3.14")

install.packages("R.utils")
remotes::install_version("tidyverse", version = "1.3.1")
remotes::install_version("Seurat", version = "3.2.0")
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github('satijalab/seurat-wrappers')
remotes::install_github("velocyto-team/velocyto.R")
remotes::install_github('cole-trapnell-lab/monocle3')
remotes::install_version("targets", version = "0.10.0")
remotes::install_version("tarchetypes", version = "0.4.1")
remotes::install_version("rmarkdown", version = "2.11")
remotes::install_version("markdown", version = "1.1")
remotes::install_version("future.callr", version = "0.7.0")
remotes::install_version("visNetwork", version = "2.1.0")
remotes::install_version("kableExtra", version = "1.3.4")
remotes::install_version("spatstat", version = "1.64-1")
remotes::install_version("Matrix", version = "1.3-3")

