# Libraries ----
library(Seurat)
library(tidyverse)
# Load data ----
seurat <- readRDS("clustered_seurat_04.rds")
# Convert to sce ----
sce <- as.SingleCellExperiment(seurat)
saveRDS(sce, "clustered_sce.rds")
