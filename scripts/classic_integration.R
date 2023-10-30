# This script integrates the data not using Harmony
# Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(sctransform)
library(svglite)
library(Matrix)
library(future)
# Load data ----
b_cluster <- readRDS("b_cluster_normalized_by_sample.rds")
load("integ_features_b")
# Integration Analysis ----
b_cluster <- PrepSCTIntegration(b_cluster, anchor.features = integ_features_b)
plan("multisession", workers = 4)
anchors <- FindIntegrationAnchors(b_cluster, normalization.method = "SCT", anchor.features = integ_features_b)
b_cluster_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

b_cluster_integrated <- RunPCA(b_cluster_integrated)
b_cluster_integrated <- RunUMAP(b_cluster_integrated, reduction = "pca", dims = 1:10)
save_plot("b_cluster_integrated.svg", b_cluster_integrated)
