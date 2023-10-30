# Importing Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(svglite)
# Loading Clustered Data ----
clustered_data <- readRDS("clustered_seurat_04.rds")
# Cluster Annotation ----
clusters_annotated <- RenameIdents(clustered_data,
                                  "0" = "B naive", # B naive
                                  "1" = "TCD4+", # TCD4+ naive
                                  "2" = "B naive", # B naive
                                  "3" = "NK", # NK
                                  "4" = "B memoria", # Desconhecido
                                  "5" = "Monocito CD14+", # Monocito CD14+
                                  "6" = "TCD8+", # TCD8+
                                  "7" = "B memoria", # TCD8+
                                  "8" = "cDC1", # cDC1
                                  "9" = "Monocito CD16+", # Monocito CD16+
                                  "10" = "cDC2", # cDC2
                                  "11" = "Plasmablastos", # Plasmablastos
                                  "12" = "Plasmablastos", # Plasmablastos
                                  "13" = "pDC", # pDC
                                  "14" = "HPC", # HPC
                                  "15" = "TCD8+ naive", # TCD8+ naive
                                  "16" = "Plaquetas", # Plaquetas
                                  "17" = "B memoria", # B memoria
                                  "18" = "Monocitos Alt", # Monocitos Alt
                                  "19" = "Plasmablastos", # Plasmablastos
                                  "20" = "B memoria", # B memoria
                                  "21" = "Treg", # Treg
                                  "22" = "NK") # NK

annotate_plot <- DimPlot(clusters_annotated, reduction = "umap", label = TRUE, label.size = 2, repel = TRUE)
save_plot("annotate_cluster_final.jpeg", annotate_plot)
save_plot("annotate_cluster_final.svg", annotate_plot)
#saveRDS(clusters_annotated, "clusters_annotated.rds")
