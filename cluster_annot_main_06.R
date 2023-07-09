# Importing Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(svglite)
# Loading Clustered Data ----
clustered_data <- readRDS("clustered_seurat_06.rds")
# Cluster Annotation ----
clusters_annotated <- RenameIdents(clustered_data,
                                   "0" = "B", "1" = "NK", "2" = "B",
                                   "3" = "B", "4" = "Monocito", "5" = "B",
                                   "6" = "Desconhecido", "7" = "Monocito", "8" = "T",
                                   "9" = "B", "10" = "T", "11" = "T",
                                   "12" = "T", "13" = "B", "14" = "DC",
                                   "15" = "B", "16" = "DC", "17" = "DC",
                                   "18" = "T", "19" = "HPC", "20" = "B",
                                   "21" = "Plaqueta", "22" = "T","23" = "B",
                                   "24" = "Monocito","25" = "B", "26" = "T",
                                   "27" = "T","28" = "B","29" = "NK", "30" = "T")

#View(clusters_annotated)
annotate_plot <- DimPlot(clusters_annotated, reduction = "umap", label = TRUE, label.size = 2, repel = TRUE)
save_plot("clusters_annotated_main_06.jpeg", annotate_plot)
save_plot("clusters_annotated_main_06.svg", annotate_plot)
saveRDS(clusters_annotated, "clusters_annotated_main_06.rds")
