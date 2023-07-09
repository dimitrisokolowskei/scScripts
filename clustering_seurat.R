# Importing Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(future)
library(svglite)
# Load Integrated Seurat Object ----
integrated_seurat <- readRDS("integrated_seurat.rds")
# Clustering Analysis ----
elbow_plot <- ElbowPlot(object = integrated_seurat, ndims = 30)
#save_plot("elbow_plot_B.svg", elbow_plot)

integrated_seurat <- FindNeighbors(object = integrated_seurat, dims = 1:10, reduction = "harmony")
plan("multisession", workers = 2)
integrated_seurat <- FindClusters(object = integrated_seurat, resolution = 1.2) # Change Cluster Resolution
cluster_plot <- DimPlot(integrated_seurat, reduction = "umap", label = FALSE, label.size = 5)
saveRDS(integrated_seurat, "clustered_seurat_12.rds")
save_plot("clustered_12.svg", cluster_plot)
