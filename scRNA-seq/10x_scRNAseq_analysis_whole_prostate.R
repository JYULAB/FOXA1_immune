# 10x Genomics scRNA-seq analysis of whole prostate tumor
# Fig 2, 3, 6

# Load Packages:
library(Seurat)
library(SeuratData)
library(scCustomize)
library(ggplot2)
library(dittoSeq)
library(MAST)
library(dplyr) 

# Read in seurat object 
seuratObj = readRDS("minimal_seurat.rds")

# color palette
polychrome_pal <- scCustomize::DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]

# add another level of cell type identities - combine granular cell types into major cell type identities - Epi, Stroma, SV, Immune
Epithelial <- c("Luminal","Spink1","Pbsn","Club Epithelial","Basal","Basal Squamous","Basal Proliferating")
SV <- c("Seminal Vesicles")
Stroma <- c("Fibroblasts","Endothelial")
Immune <- c("Monocytes","Macrophage","B cells","T cells")

cl.mat <- cbind(c(rep("Epithelial", length(Epithelial)), rep("Stroma", length(Stroma)), rep("Immune", length(Immune)), rep("SV", length(SV))
),
c(Epithelial, Stroma, Immune, SV))

cl.vec <- seuratObj$celltypes_granular
ct.vec <- rep(NA, length(cl.vec))
for(x in unique(cl.mat[,1])){
  cl.x <- cl.mat[cl.mat[,1]==x,2]
  ct.vec[which(cl.vec%in%cl.x)] <- x
}
seuratObj$celltypes_simple <- ct.vec
seuratObj$celltypes_simple <- factor(seuratObj$celltypes_simple, levels = c('Epithelial', 'Stroma','Immune','SV'))
DimPlot(seuratObj, group.by = 'celltypes_simple', cols = polychrome_pal)

#change order of clusters in seurat object
seuratObj$orig.ident <- factor(seuratObj$orig.ident, levels = c("P12","PF12","P18","PF18","XP","XPF"))

# subset on intact mice 
intact_seuratObj <- seuratObj[,seuratObj$timepoint != "X"]

# Fig 2A
# cell type umaps - dimplot()
DimPlot(intact_seuratObj, group.by = "orig.ident", cols = polychrome_pal, shuffle = T)
DimPlot(intact_seuratObj, group.by = "celltypes_coarse", cols = polychrome_pal, shuffle = T)
DimPlot(intact_seuratObj, group.by = "celltypes_granular", cols = polychrome_pal, shuffle = T)


#-------------------------------------------------------------------------------
# Determining cell type proportions using DittoSeq

# Remove SV for cell type proportions
seuratObj_noSV <- intact_seuratObj[,intact_seuratObj$celltypes_granular != "Seminal Vesicles"]

# all cell types - reordered
dittoBarPlot(
  object = seuratObj_noSV,
  var = "celltypes_coarse",
  group.by = "orig.ident",
  color.panel = polychrome_pal,
  var.labels.reorder = c(3,2,6,4,5,7,1),
  x.reorder = c(1,3,2,4))

#-------------------------------------------------------------------------------
# Gene lists
major_markers_granular <- c("Epcam","Krt8","Wfdc2","Spink1","Pbsn","Tff3","Pate4","Krt5","Krt14","Trp63","Krt6a", "Mki67","Lum","Ackr1","S100a8","Lyz2","Ms4a1","Cd3e","Nkg7")
major_markers_coarse <- c("Epcam","Krt8","Krt5","Trp63","Krt6a","Pate4","Lum","S100a8","Lyz2","Ms4a1","Cd3e","Nkg7")

# Violin plots of major cell type markers - Fig 2
# coarse cell types
VlnPlot(intact_seuratObj, features = major_markers_coarse, group.by = "celltypes_coarse", stack = TRUE,  flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Cell type markers")

# granular cell types
VlnPlot(intact_seuratObj, features = major_markers_granular, group.by = "celltypes_granular", stack = TRUE,  flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Cell type markers")

#-------------------------------------------------------------------------------
# Dot plot of gene expression within epithelial cells 

# subset on Epithelial
Epi_seuratObj <- seuratObj_noSV[,seuratObj_noSV$celltypes_simple == "Epithelial"]
DimPlot(Epi_seuratObj, group.by = "celltypes_simple", cols = polychrome_pal)

# gene lists
basal_squam <- c("Trp63", "Krt5","Krt14","Lgals7", "Krt6a","Krt13","Krt16")
luminal <- c("Foxa1","Ar","Krt8","Krt18","Cd24a")

# Fig 2
Idents(Epi_seuratObj) <- Epi_seuratObj$orig.ident
DotPlot_scCustom(Epi_seuratObj, features = basal_squam, group.by = "orig.ident") + RotatedAxis() + labs(title=" Basal-squamous")
DotPlot_scCustom(Epi_seuratObj, features = luminal, group.by = "orig.ident") + RotatedAxis() + labs(title="     Luminal")
DotPlot_scCustom(Epi_seuratObj, features = c("Vim","Snai2","Twist1"), group.by = "orig.ident") + RotatedAxis() + labs(title="     EMT") + scale_size(range = c(1.7, 6))

# Fig 6
DotPlot_scCustom(Epi_seuratObj, features = c("Tgfb3","Tgfb2","Tgfb1"), group.by = "orig.ident") + RotatedAxis() + labs(title="Epithelial")

#-------------------------------------------------------------------------------
# Feature plots

# Fig 2
FeaturePlot_scCustom(seurat_object = intact_seuratObj, features = "Krt8", split.by = "orig.ident", order = T, pt.size = 0.2)
FeaturePlot_scCustom(seurat_object = intact_seuratObj, features = "Krt5", split.by = "orig.ident", order = T, pt.size = 0.2)
FeaturePlot_scCustom(seurat_object = intact_seuratObj, features = "Krt6a", split.by = "orig.ident", order = T, pt.size = 0.2)
FeaturePlot_scCustom(seurat_object = intact_seuratObj, features = major_markers_granular)

# Fig 6
FeaturePlot_scCustom(seurat_object = intact_seuratObj, features = c("Tgfb3"), split.by = "orig.ident", order = T, pt.size = 0.2)
FeaturePlot_scCustom(seurat_object = intact_seuratObj, features = c("Tgfb2"), split.by = "orig.ident", order = T, pt.size = 0.2)
FeaturePlot_scCustom(seurat_object = intact_seuratObj, features = c("Tgfb1"), split.by = "orig.ident", order = T, pt.size = 0.2)
FeaturePlot_scCustom(seurat_object = intact_seuratObj, features = c("Tgfbr2"), split.by = "orig.ident", order = T, pt.size = 0.2)
FeaturePlot_scCustom(seurat_object = intact_seuratObj, features = c("Tgfbr1"), split.by = "orig.ident", order = T, pt.size = 0.2)

#-------------------------------------------------------------------------------
# Differential gene expression in PF vs P epithelial cells - using MAST
# Fig 3A-B

Idents(Epi_seuratObj) <- Epi_seuratObj$orig.ident
w18_Epi_DEG_MAST <- FindMarkers(
  Epi_seuratObj,
  ident.1 = "PF18",
  ident.2 = "P18", 
  min.pct = 0.15, logfc.threshold = -Inf, test.use = "MAST")
write.table(w18_Epi_DEG_MAST, file = "w18_Epi_DEG_MAST.csv", sep = ",", quote = FALSE)

Idents(Epi_seuratObj) <- Epi_seuratObj$orig.ident
w12_Epi_DEG_MAST <- FindMarkers(
  Epi_seuratObj,
  ident.1 = "PF12",
  ident.2 = "P12", 
  min.pct = 0.15, logfc.threshold = -Inf, test.use = "MAST")
write.table(w12_Epi_DEG_MAST, file = "w12_Epi_DEG_MAST.csv", sep = ",", quote = FALSE)

