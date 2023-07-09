# Libraries ----
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
# Load data ----
sce <- readRDS("cluster_annotated_singler.rds")
# Convert to sce ----
seurat <- as.Seurat(sce)
saveRDS(seurat, "clustered_sce_to_seurat.rds")
