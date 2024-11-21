# Fig 6 - CellChat analysis of 10x scRNA-seq whole prostate tumor data

# citation
#  title: "Comparison analysis of multiple datasets using CellChat"
# author: "Suoqin Jin"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#  html_document:
#  toc: true
#theme: united
#mainfont: Arial
#vignette: >
#  %\VignetteIndexEntry{Comparison analysis of multiple datasets using CellChat}
#%\VignetteEngine{knitr::rmarkdown}
#%\VignetteEncoding{UTF-8}

library(CellChat)
library(patchwork)
library(Seurat)
library(draw)
library(ComplexHeatmap)
library(circlize)
library(export)
library(dplyr) 
library(tidyverse)
options(stringsAsFactors = FALSE)

# Create CellChat object
# read in seurat object
seuratObj = readRDS("minimal_seurat.rds")

data.input = GetAssayData(seuratObj, slot = "data", assay = "RNA")
meta = seuratObj@meta.data
meta$labels = meta$celltypes_coarse
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
unique(meta$orig.ident)
unique(meta$labels)

#P18
cell.use_P18 = rownames(meta)[meta$orig.ident == "P18"]
data.input_P18 = data.input[, cell.use_P18]
meta_P18 = meta[cell.use_P18,]
unique(meta_P18$labels)
cellchat_P18 <- createCellChat(object = data.input_P18, meta = meta_P18, group.by = "labels")
cellchat_P18 <- updateCellChat(cellchat_P18)

#PF18
cell.use_PF18 = rownames(meta)[meta$orig.ident == "PF18"]
data.input_PF18 = data.input[, cell.use_PF18]
meta_PF18 = meta[cell.use_PF18,]
unique(meta_PF18$labels)
cellchat_PF18 <- createCellChat(object = data.input_PF18, meta = meta_PF18, group.by = "labels")
cellchat_PF18 <- updateCellChat(cellchat_PF18)

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis - secreted signaling
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 

# set the used database for each object
cellchat_P18@DB <- CellChatDB.use
cellchat_PF18@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes to save computation cost
cellchat_P18<- subsetData(cellchat_P18) 
cellchat_PF18 <- subsetData(cellchat_PF18) 

cellchat_P18<- identifyOverExpressedGenes(cellchat_P18)
cellchat_P18<- identifyOverExpressedInteractions(cellchat_P18)
cellchat_PF18 <- identifyOverExpressedGenes(cellchat_PF18)
cellchat_PF18 <- identifyOverExpressedInteractions(cellchat_PF18)

#-------------------------------------------------------------------------------
#Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network
cellchat_P18<- computeCommunProb(cellchat_P18, population.size = TRUE) #default 25% truncated mean
cellchat_PF18 <- computeCommunProb(cellchat_PF18, population.size = TRUE) #default 25% truncated mean

df.net_P18<- subsetCommunication(cellchat_P18)
df.net_PF18 <- subsetCommunication(cellchat_PF18)

cellchat_P18<- computeCommunProbPathway(cellchat_P18)
cellchat_PF18 <- computeCommunProbPathway(cellchat_PF18)

cellchat_P18<- aggregateNet(cellchat_P18)
cellchat_PF18 <- aggregateNet(cellchat_PF18)

cellchat_P18<- updateCellChat(cellchat_P18)
cellchat_PF18 <- updateCellChat(cellchat_PF18)

#-------------------------------------------------------------------------------
# Visualization of cell-cell communication network

pathways.show.tgfb <- c("TGFb") 
pathways.show.MIF <- c("MIF") 
pathways.show.PTN <- c("PTN") 

# circle plot
netVisual_aggregate(cellchat_PF18, signaling = pathways.show.tgfb, layout = "circle")
netVisual_aggregate(cellchat_PF18, signaling = pathways.show.PTN, layout = "circle")
netVisual_aggregate(cellchat_PF18, signaling = pathways.show.MIF, layout = "circle")

# Fig 6C
# Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)

# outgoing patterns
selectK(cellchat_PF18, pattern = "outgoing")
nPatterns = 5
cellchat_PF18 <- identifyCommunicationPatterns(cellchat_PF18, pattern = "outgoing", k = nPatterns)
netAnalysis_dot(cellchat_PF18, pattern = "outgoing")

#-------------------------------------------------------------------------------
# COMPARING MULTIPLE CELLCHAT OBJECTS
# merge cellchat objects
object.list <- list(P18 = cellchat_P18, PF18 = cellchat_PF18)
cellchat_w18 <- mergeCellChat(object.list, add.names = names(object.list))

# Predict general principles of cell-cell communication
# compare overall # and strength of interactions
gg1 <- compareInteractions(cellchat_w18, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_w18, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# Differential number of interactions among major cell types
group.cellType <- c(rep("Epi", 3), rep("Stroma", 1), rep("Immune", 4))
group.cellType <- factor(group.cellType, levels = c("Epi", "Stroma", "Immune"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat_w18  <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_diffInteraction(cellchat_w18, weight.scale = T, measure = "count.merged", label.edge = T) # Fig 6A

# Identify signaling pathways enriched in PF vs P
rankNet(cellchat_w18, mode = "comparison", stacked = T, do.stat = TRUE) 
rankNet(cellchat_w18, mode = "comparison", stacked = F, do.stat = TRUE) # Fig 6B