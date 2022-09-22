install.packages("https://cran.r-project.org/src/contrib/remotes_2.4.2.tar.gz")
remotes::install_github("cole-trapnell-lab/leidenbase", upgrade = "never")
remotes::install_version("spatstat", version = "1.64-1")
remotes::install_github("satijalab/seurat", ref = "84e5d20475f4d7db6b5b7d08ed5d4317831f662a", upgrade = "never")

install.packages("BiocManager")
BiocManager::install(version = '3.12', update = FALSE, ask = FALSE)
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
  version = "3.12", update = FALSE)

install.packages("R.utils")
remotes::install_version("tidyverse", version = "1.3.0")
remotes::install_github("mojaveazure/loomR", upgrade = "never")
remotes::install_github('cole-trapnell-lab/monocle3', upgrade = "never", ref = "1.0.0")
remotes::install_version("targets", version = "0.1.0")
remotes::install_version("tarchetypes", version = "0.1.0")
remotes::install_version("rmarkdown", version = "2.7")
remotes::install_version("markdown", version = "1.1")
remotes::install_version("visNetwork", version = "2.0.9")
remotes::install_version("kableExtra", version = "1.3.4")
remotes::install_version("sctransform", version = "0.3.2")
