# Libraries ----
library(SingleR)
library(SingleCellExperiment)
# Load Reference and query data ----
reference <- readRDS("reference_mid.rds") # 
sce <- readRDS("clustered_sce.rds") # Already converted Seurat --> sce 
# SingleR analysis ----
pred <- SingleR(test = sce, ref = reference, 
                labels = reference$label.fine)

sce$singler <- pred$labels
saveRDS(sce, "cluster_annotated_singler_fine.rds")

