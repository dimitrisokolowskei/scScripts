# libraries ----
library(Seurat)
library(tidyverse)
library(CellChat)
library(patchwork)
library(future)
library(NMF)
library(ggalluvial)
library(cowplot)
library(svglite)
# import data ----
time <- readRDS("time7.rds")
# Prepare cellchat objects ----
data <- GetAssayData(time, assay = "SCT", slot = "data") # normalized data matrix
labels <- Idents(time)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
# Create cellchat object ----
cellchat <- createCellChat(object = data)
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
# Preprocessing for CCC analysis ----
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
levels(cellchat@idents)
#View(cellchat)
future::plan("multisession", workers = 4)
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# Communication Probability Computation ----
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# Visualization ----
groupSize <- as.numeric(table(cellchat@idents))
net_interactions <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
net_strenghts <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- cellchat@net$weight
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Visualize Sgnalling Pathway ----
par(mfrow=c(1,1))
cellchat@netP$pathways # See pathways
pathways.show <- c("GALECTIN") 
vertex.receiver = seq(1,4) # a numeric vector. 
net_aggreg_hire <- netVisual_aggregate(cellchat, signaling = "CXCL", layout = "hierarchy", vertex.receiver = vertex.receiver, vertex.size = 2)
# Circle plot
#par(mfrow=c(1,1))
#netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
pathways.show <- c("BAG") 
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# Gene Expression 
plotGeneExpression(cellchat, signaling = "BAG")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair ----
net_contribution <- netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # Modificar 
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
net_individual <- netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")


# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways ----
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 16, targets.use = c(1:16), remove.isolate = FALSE) # Change sources 

# Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling ----
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
cellchat@netP$pathways # See pathways
pathways.show <- c("CD30") 
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 10, height = 3, font.size = 10)
signal_role_
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

signal_role_tgfb
# Visualize the dominant senders (sources) and receivers (targets) in a 2D space ----
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together ----
# Outgoing patterns
selectK(cellchat, pattern = "outgoing") 
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# Ingoing Patterns
selectK(cellchat, pattern = "incoming") 
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
