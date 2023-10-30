# This script splits seurat Object samples
# Importing Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(sctransform)
# Load filtered Seurat Object
b_cluster <- readRDS("b_cluster_subset_by_sample.rds") # seurat_filtered_no_doublet_by_score.RData

# Splitting data sample ----
split_seurat_sample <- SplitObject(b_cluster, split.by = "sample")

# Normalizing "samples" individualy ----
for (i in 1:length(split_seurat_sample)) {
  split_seurat_sample[[i]] <- SCTransform(split_seurat_sample[[i]], vars.to.regress = c("mitoRatio"))
} # Normalize without mitoration regression

saveRDS(split_seurat_sample, file="b_cluster_subset_normalized_by_sample.rds")
