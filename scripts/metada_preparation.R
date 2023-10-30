# Brief Description ----
# I perform quality control to remove cell with undesired features and removing doublets.
# Importing Libraries ----
library(SeuratObject)
library(Seurat)
library(tidyverse)
# Import Merged Data ----
load("merged_seurat.RData")
# Analysis
#View(merged_seurat@meta.data)
# Compute UMI and Mito Ratio ----
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
# Add aditional columns to metadata
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
# Add sample ID
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^2047_"))] <- "2047"
metadata$sample[which(str_detect(metadata$cells, "^2049_"))] <- "2049"
metadata$sample[which(str_detect(metadata$cells, "^2051_"))] <- "2051"
metadata$sample[which(str_detect(metadata$cells, "^2052_"))] <- "2052"
metadata$sample[which(str_detect(metadata$cells, "^2053_"))] <- "2053"
metadata$sample[which(str_detect(metadata$cells, "^2055_"))] <- "2055"
metadata$sample[which(str_detect(metadata$cells, "^P1_"))] <- "P1"
metadata$sample[which(str_detect(metadata$cells, "^P2_"))] <- "P2"
metadata$sample[which(str_detect(metadata$cells, "^P3_"))] <- "P3"
metadata$sample[which(str_detect(metadata$cells, "^P4_"))] <- "P4"
metadata$sample[which(str_detect(metadata$cells, "^P5_"))] <- "P5"
metadata$sample[which(str_detect(metadata$cells, "^P6_"))] <- "P6"
metadata$sample[which(str_detect(metadata$cells, "^P7_"))] <- "P7"
metadata$sample[which(str_detect(metadata$cells, "^P8_"))] <- "P8"
metadata$sample[which(str_detect(metadata$cells, "^P9_"))] <- "P9"
#View(metadata)
# Add Time periods
metadata$time <- NA
metadata$time[which(str_detect(metadata$cells, "^2047_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^2049_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^2051_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^2052_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^2053_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^2055_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^P1_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^P2_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^P3_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^P4_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^P5_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^P6_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^P7_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^P8_GEX0"))] <- "T0"
metadata$time[which(str_detect(metadata$cells, "^P9_GEX0"))] <- "T0"

metadata$time[which(str_detect(metadata$cells, "^2047_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^2049_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^2051_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^2052_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^2053_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^2055_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^P1_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^P2_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^P3_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^P4_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^P5_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^P6_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^P7_GEX7"))] <- "T7"
metadata$time[which(str_detect(metadata$cells, "^P8_GEX7"))] <- "T7"
#metadata$time[which(str_detect(metadata$cells, "^P9_GEX7"))] <- "T7"

metadata$time[which(str_detect(metadata$cells, "^2047_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^2049_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^2051_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^2052_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^2053_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^2055_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^P1_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^P2_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^P3_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^P4_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^P5_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^P6_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^P7_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^P8_GEX21"))] <- "T21"
metadata$time[which(str_detect(metadata$cells, "^P9_GEX21"))] <- "T21"

metadata$time[which(str_detect(metadata$cells, "^2047_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^2049_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^2051_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^2052_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^2053_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^2055_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^P1_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^P2_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^P3_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^P4_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^P5_GEX28"))] <- "T28"
#metadata$time[which(str_detect(metadata$cells, "^P6_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^P7_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^P8_GEX28"))] <- "T28"
metadata$time[which(str_detect(metadata$cells, "^P9_GEX28"))] <- "T28"
# Add study ID
metadata$study <- NA
metadata$study[which(str_detect(metadata$cells, "^20"))] <- "Arunachalam et al., 2021"
metadata$study[which(str_detect(metadata$cells, "^P"))] <- "Brewer et al., 2021"
# Rename Columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Assign metadata back and generate .RData file
merged_seurat@meta.data <- metadata
save(merged_seurat, file="merged_metadata_seurat.RData") # .RData with metadata prepared

#View(metadata)
