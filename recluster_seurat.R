# Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(svglite)
# Data loading ----
clustered_data <- readRDS("clusters_annotated_main_only.rds") # Standard File: clusters_annotated.rds
# Recluster B cells ----
b_cluster <- subset(clustered_data, idents = "T", invert = FALSE)
#View(b_cluster)
saveRDS(b_cluster, "T_cluster_from_main_only.rds")
