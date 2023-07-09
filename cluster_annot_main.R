# This script performs UMAP with only the main cell types
# Importing Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(svglite)
# Loading Clustered Data ----
clustered_data <- readRDS("clustered_seurat_04.rds")
# Cluster Annotation ----
clusters_annotated <- RenameIdents(clustered_data,
                                  "0" = "B", # B naive
                                  "1" = "T", # TCD4+ naive
                                  "2" = "B", # B naive
                                  "3" = "NK", # NK
                                  "4" = "B", # Desconhecido
                                  "5" = "Monocito", # Monocito CD14+
                                  "6" = "T", # TCD8+
                                  "7" = "B", # B memoria
                                  "8" = "DC", # cDC1
                                  "9" = "Monocito", # Monocito CD16+
                                  "10" = "DC", # cDC2
                                  "11" = "Plasmablastos", # Plasmablastos
                                  "12" = "Plasmablastos", # Plasmablastos
                                  "13" = "DC", # pDC
                                  "14" = "HPC", # HPC
                                  "15" = "T", # TCD8+ naive
                                  "16" = "Plaquetas", # Plaquetas
                                  "17" = "B", # B memoria
                                  "18" = "Monocito", # Monocitos Alt
                                  "19" = "Plasmablastos", # Plasmablastos
                                  "20" = "B", # B memoria
                                  "21" = "T", # Treg
                                  "22" = "NK") # NK

annotate_plot <- DimPlot(clusters_annotated, reduction = "umap", label = TRUE, label.size = 2, repel = TRUE)
save_plot("annotate_cluster_main.svg", annotate_plot)
save_plot("annotate_cluster_main.jpeg", annotate_plot)
saveRDS(clusters_annotated, "clusters_annotated_main_only.rds")
