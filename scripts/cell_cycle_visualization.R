# This script normalizes and Visualize cell cycle differences 
# Importing Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(svglite)
# Load filtered Seurat Object
cluster <- readRDS("clusters_annotated.rds")
# Regression Analysis (Cell cycle and Mito) ----
exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, 
                      as.is = TRUE, row.names = 1) # Importing file with cell cycle genes

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat_phase <- CellCycleScoring(cluster, 
                                 g2m.features = g2m.genes, 
                                 s.features = s.genes)

# Plot the PCA colored by cell cycle phase
pca_img <- DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase",
        raster =FALSE) 

saveRDS(seurat_phase, file = "clustered_seurat_cellcycle_annotation.rds")
save_plot("cell_cycle.jpeg", pca_img)
save_plot("cell_cycle.svg", pca_img)
