# Importing Library
library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(svglite)
# Load Data imputated ----
merged_seurat_again <- readRDS("merged_b_cluster_seurat_again.rds") # Data normalized and merged again (Ex: merged_{optional}_again.rds)
load("integ_features_b")
# Integration Analysis ----
DefaultAssay(merged_seurat_again) <- "SCT"

VariableFeatures(merged_seurat_again) <- integ_features_b

merged_seurat <- RunPCA(merged_seurat_again, assay = "SCT", npcs = 30) # Assay SCT or arla
#non_integ_umap <- RunUMAP(merged_seurat, assay = "SCT", dims = 1:30)
#umap_plot_noninteg <- DimPlot(non_integ_umap, group.by="study")
#save_plot("umap_plot_noninteg.svg", umap_plot_noninteg)
#pca_plot <- PCAPlot(merged_seurat, split.by = "sample")
#save_plot("PCA.jpeg", pca_plot)

harmonized_seurat <- RunHarmony(merged_seurat, 
                                group.by.vars = c("sample", "study"), 
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:30)
saveRDS(harmonized_seurat, file="integrated_B_seurat.rds")

umap_plot_sample <- DimPlot(harmonized_seurat, group.by="sample")
#umap_plot_study <- DimPlot(harmonized_seurat, group.by="study")
#umap_plot_time <- DimPlot(harmonized_seurat, group.by="time")
save_plot("UMAP_B_cluster_subset.svg", umap_plot_sample)
#save_plot("UMAP_study.svg", umap_plot_study)
#save_plot("UMAP_time.svg", umap_plot_time)

