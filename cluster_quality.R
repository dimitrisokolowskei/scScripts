# Importing Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(svglite)
# Load Integrated Seurat Object ----
clusters <- readRDS("clustered_seurat_cellcycle_annotation.rds") # Import Cluster with specific resolution
# Cluster Quality plotting ----
n_cells <- FetchData(clusters, 
                     vars = c("ident", "orig.ident")) %>%
        dplyr::count(ident, orig.ident) %>%
        tidyr::spread(ident, n)


# UMAP of cells in each condition (sample, time or study) ----
umap_by_sample <- DimPlot(clusters, label = TRUE, split.by = "study")  + NoLegend()
save_plot("umap_by_study.svg", umap_by_sample)
save_plot("umap_by_study.jpeg", umap_by_sample)

# Explore whether clusters segregate by cell cycle phase
umap_by_cycle <- DimPlot(clusters, label = TRUE, split.by = "Phase")  + NoLegend()
save_plot("umap_by_cycle.svg", umap_by_cycle)
save_plot("umap_by_cycle.jpeg", umap_by_cycle)

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
clusters_qc_features <- FeaturePlot(clusters, reduction = "umap", 
                                    features = metrics, pt.size = 0.4, 
                                    order = TRUE, min.cutoff = 'q10',
                                    label = TRUE)
save_plot("clusters_qc_features.svg", clusters_qc_features)
save_plot("clusters_qc_features.jpeg", clusters_qc_features)
