# This script identified and annotated doublets in Seurat Object metadata
# Import Libraries ----
library(Seurat)
library(SeuratObject)
library(scDblFinder)
library(Matrix)
library(BiocParallel)
# Load Data ----
load("merged_metadata_seurat.RData")
# Doublet Identification ----
set.seed(123)
sce <- as.SingleCellExperiment(merged_seurat)
bp <- MulticoreParam(2, RNGseed=1234)
sce <- scDblFinder(sce, samples="sample", BPPARAM=bp)
seurat_filtered <- as.Seurat(sce)
save(seurat_filtered, file="merged_metadata_doublet_annotated_seurat.RData")
