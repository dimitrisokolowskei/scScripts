# Import Libraries ----
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(Matrix)
# Load Data ----
load("merged_metadata_doublet_annotated_seurat.RData")
#View(seurat_filtered@meta.data)
#dim(seurat_filtered)

before90 <- table(seurat_filtered@meta.data$scDblFinder.score >= 0.80)
# Quality Filtering Analysis ----
filtered_seurat <- subset(x = seurat_filtered, 
                          subset= (nUMI >= 600) & 
                            (nGene >= 350) &
                            (nGene <= 9000) &
                            (log10GenesPerUMI > 0.82) & 
                            (mitoRatio < 0.10) &
                            (scDblFinder.score < 0.80)) 


#dim(filtered_seurat)
#View(filtered_seurat@meta.data)

counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
#dim(filtered_seurat)
save(filtered_seurat, file="seurat_filtered_no_doublet_by_score.RData")
