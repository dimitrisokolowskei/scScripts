# Importing Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(svglite)
# Loading Clustered Data ----
clustered_data <- readRDS("clustered_seurat_06.rds")
# Cluster Annotation ----
clusters_annotated <- RenameIdents(clustered_data,
                                   "0" = "B naive", "1" = "NK", "2" = "B naive", 
                                   "3" = "B memoria", "4" = "Monocito CD14+", "5" = "B memoria", 
                                   "6" = "Desconhecido", "7" = "Monocito CD16+", "8" = "TCD4+", 
                                   "9" = "B naive", "10" = "TCD8+", "11" = "TCD8+", 
                                   "12" = "Treg", "13" = "B memoria", "14" = "cDC2", 
                                   "15" = "Plasmablastos", "16" = "cDC1", "17" = "pDC", 
                                   "18" = "TCD4+", "19" = "HPC", "20" = "B memoria", 
                                   "21" = "Plaqueta", "22" = "TCD4+","23" = "B memoria", 
                                   "24" = "Monocito Alt","25" = "Plasmablastos", "26" = "TCD4+", 
                                   "27" = "TCD4+","28" = "B memoria","29" = "NK", "30" = "TCD4+") 

View(clusters_annotated)
annotate_plot <- DimPlot(clusters_annotated, reduction = "umap", label = TRUE, label.size = 2, repel = TRUE)
save_plot("clusters_annotated_06.jpeg", annotate_plot)
save_plot("clusters_annotated_06.svg", annotate_plot)
saveRDS(clusters_annotated, "clusters_annotated_06.rds")