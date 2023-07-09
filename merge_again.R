# Script for performing data merging after individual SCTransform samples normalization
# Importing Library
library(Seurat)
library(tidyverse)
library(cowplot)
# Load splited and Normalized merged data ----
split_seurat_sample <- readRDS("b_cluster_subset_normalized_by_sample.rds")
#View(split_seurat_sample)
# Analysis ----
integ_features_b_subset <- SelectIntegrationFeatures(object.list = split_seurat_sample, nfeatures = 3000)

merged_seurat <- merge(x = split_seurat_sample[[1]],
                       y = split_seurat_sample[2:length(split_seurat_sample)],
                       merge.data = TRUE)

#View(merged_seurat)

saveRDS(merged_seurat, file="merged_b_cluster_subset_seurat_again.rds")
save(integ_features_b_subset, file="integ_features_b_subset")
