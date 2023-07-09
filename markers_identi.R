# Importing Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(metap)
library(multtest)
# Loading Clusterized Data ----
cluster_data <- readRDS("clustered_seurat_06.rds")
# Markers Identification ----
DefaultAssay(cluster_data) <- "SCT"

#markers <- FindConservedMarkers(cluster_data,
#                                ident.1 = 1,
#                                grouping.var = "sample",
#                                only.pos = TRUE,
#                                logfc.threshold = 1)

#write.csv(markers, "markers_logfc2.csv")

# Find all markers in each cluster ----
sct_prepared <- PrepSCTFindMarkers(cluster_data)



markers_all <- FindAllMarkers(sct_prepared, 
                              only.pos = TRUE,
                              logfc.threshold = 0.10)   

#View(markers)
write.csv(markers_all, "markers_clustered_06.csv")



