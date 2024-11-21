# PIP-seq scRNA-seq analysis of CD45+ live cells
# Fig 4

# Load Packages:
library(Seurat)
library(SeuratData)
library(scCustomize)
library(ggplot2)
library(dittoSeq)
library(MAST)
library(dplyr) 

# Read in seurat object
pip_seuratObj = readRDS("pip-seq_integrated.RDS")

#color palette
polychrome_pal <- scCustomize::DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]

# add another level of cell type identities - combine macrophage subsets into one general group - cell_types_with_t_coarse_v2
Macrophage <- c("M2 Macrophage","Macrophage")
B_cell <- c("B Cell")
CD8_T <- c("CD8 T Cell")
DC <- c("Dendritic Cell")
Monocyte <- c("Monocyte")
NK <- c("NK")
T_Helper <- c("T helper")
Tumor_Epi <- c("Tumor Epithelial")

cl.mat <- cbind(c(rep("Macrophage", length(Macrophage)), rep("B Cell", length(B_cell)), rep("CD8 T", length(CD8_T)), 
                rep("DC", length(DC)), rep("Monocyte", length(Monocyte)) , rep("NK", length(NK)), 
                rep("T Helper", length(T_Helper)), rep("Tumor Epithelial", length(Tumor_Epi))),
c(Macrophage, B_cell, CD8_T,DC, Monocyte, NK, T_Helper, Tumor_Epi ))

cl.vec <- pip_seuratObj$cell_types_with_t_coarse
ct.vec <- rep(NA, length(cl.vec))
for(x in unique(cl.mat[,1])){
  cl.x <- cl.mat[cl.mat[,1]==x,2]
  ct.vec[which(cl.vec%in%cl.x)] <- x
}
pip_seuratObj$cell_types_with_t_coarse_v2 <- ct.vec
pip_seuratObj$cell_types_with_t_coarse_v2 <- factor(pip_seuratObj$cell_types_with_t_coarse_v2, 
                                                    levels = c('Macrophage', 'B Cell', 'CD8 T','DC', 'Monocyte', 'NK', 'T Helper', 'Tumor Epithelial'))
DimPlot(pip_seuratObj, group.by = 'cell_types_with_t_coarse_v2', cols = polychrome_pal, shuffle = T)

#change order of clusters in seurat object
pip_seuratObj$cell_types_with_t_coarse_v2 <- factor(pip_seuratObj$cell_types_with_t_coarse_v2, levels = c("Monocyte", "Macrophage","DC", "B Cell","CD8 T","T Helper","NK","Tumor Epithelial"))
DimPlot(pip_seuratObj, group.by = 'cell_types_with_t_coarse_v2', cols = polychrome_pal, shuffle = T)

#change order of clusters in seurat object
pip_seuratObj$cell_types_with_t_granular <- factor(pip_seuratObj$cell_types_with_t_granular, levels = c("Monocyte", "Macrophage","M2 Macrophage","Dendritic Cell", "B Cell","Tregs","Naive/Memory T helper","Th17","Cycling CD8","CTL3","CTL2","CTL1","NK","Tumor Epithelial"))
DimPlot(pip_seuratObj, group.by = 'cell_types_with_t_granular', cols = polychrome_pal, shuffle = T)

# remove contaminating tumor epithelial cells
immune_seuratObj <- pip_seuratObj[,pip_seuratObj$cell_types_with_t_coarse_v2 != "Tumor Epithelial"]
DimPlot(immune_seuratObj, group.by = 'cell_types_with_t_coarse_v2', cols = polychrome_pal, shuffle = T)

#subset on macrophages
mac_seuratObj <- pip_seuratObj[,pip_seuratObj$cell_types_with_t_coarse_v2 == "Macrophage"]

#subset T cells
Tcell_seuratObj <- pip_seuratObj[,pip_seuratObj$cell_types_with_t_coarse_v2 == c("T Helper","CD8 T")]

#subset on CD8 T cells
CD8_seuratObj <- pip_seuratObj[,pip_seuratObj$cell_types_with_t_coarse_v2 == "CD8 T"]

# Fig 4A - umap visualizations 
DimPlot(immune_seuratObj, group.by = 'cell_types_with_t_coarse_v2', cols = polychrome_pal, shuffle = T)
DimPlot(immune_seuratObj, group.by = 'orig.ident', cols = polychrome_pal, shuffle = T)

#-------------------------------------------------------------------------------
# Cell type markers

# Gene lists 
coarse_immune_markers <- c("Itgam","S100a8",	"S100a9","Cd14","Adgre1",	"Mafb","Lyz2","C1qb", "Itgax", "Flt3", "Cd209a","H2-DMb2"	,	"Ms4a1", "Cd19","Cd3e",	"Cd8a",	"Cd8b1","Gzmb","Cd4",	"Foxp3", "Pdcd1","Ctla4","Klrd1","Nkg7")
Tcell_granular_markers <- c("Cd8a",	"Cd8b1","Gzma",	"Gzmb",	"Gzmm",	"Gzmk", "Prf1", "Havcr2",	"Tigit","Tox","Ifng","Tnf", "Top2a","Cdk1","Mki67","Cd4","Rorc",	"Rora",	"Il17a",	"Il7r", "Sell",	"Ccr7",	"Lef1",	"Il4ra","Foxp3", "Pdcd1","Ctla4")
CytotoxicityScore <-	c("Gzma",	"Gzmb",	"Gzmm",	"Gzmk",	"Prf1")
M2_score <- c("Msr1",	"Mrc1",	"Cd163",  "Il10",	"Ccl17",	"Ccl22",	"Ccl24",	"Arg1",	"Cd200r1",	"Pdcd1lg2",	"Cd274",	"Csf1r",	"Il1rn",	"Il1r2",	"Il4ra",	"Ccl20",	"Lyve1",	"Vegfa",	"Vegfb",	"Vegfc",	"Vegfd",	"Egf",	"Tgfb1",	"Tgfb2",	"Tgfb3",	"Mmp14",	"Mmp19",	"Mmp9",	"Wnt7b",	"Fasl",	"Tnfsf12",	"Tnfsf8",	"Cd276",	"Vtcn1",	"Fn1",	"Irf4")
M1_score <- c("Nos2",	"Il12a",	"Fcgr1","Cd80",	"Il23a",	"Cxcl9",	"Cxcl10",		"Cd86",	"Il1a",	"Il1b",	"Il6",	"Tnf",	"Cd68",	"Ccl5",	"Irf5",	"Irf1",	"Cd40",	"Ido1",	"Kynu",	"Ccr7")

# dot plots of cell type markers 
Idents(immune_seuratObj) <- immune_seuratObj$cell_types_with_t_coarse_v2
DotPlot_scCustom(immune_seuratObj, features = coarse_immune_markers, group.by = "cell_types_with_t_coarse_v2") + RotatedAxis()+ scale_colour_gradient2(low = "blue", mid = "white", high = "red") + scale_size(range = c(1,6))

Idents(Tcell_seuratObj) <- Tcell_seuratObj$cell_types_with_t_granular
DotPlot_scCustom(Tcell_seuratObj, features = Tcell_granular_markers, group.by = "cell_types_with_t_granular")+ scale_size(range = c(1,6))   + scale_colour_gradient2(low = "blue", mid = "white", high = "red") + RotatedAxis()

# feature plots of major cell type markers
FeaturePlot_scCustom(seurat_object = immune_seuratObj, features = c("Itgax","Adgre1","S100a9","Cd209a","Cd3e","Cd8a","Cd4", "Nkg7","Ms4a1","Rorc"))

#-------------------------------------------------------------------------------
# Macrophages
# Add module score - M1 M2
mac_seuratObj <- AddModuleScore(mac_seuratObj,
                                features = list(M1_score),
                                name="M1_score")
VlnPlot(mac_seuratObj, features = "M1_score1", group.by = "orig.ident") + stat_summary(fun.y = median, geom='point', size = 15, colour = "white", shape = 95)

mac_seuratObj <- AddModuleScore(mac_seuratObj,
                                features = list(M2_score),
                                name="M2_score")
VlnPlot(mac_seuratObj, features = "M2_score1", group.by = "orig.ident") +  stat_summary(fun.y = median, geom='point', size = 15, colour = "white", shape = 95) 

#-------------------------------------------------------------------------------
# T cells
# Dimplot
DimPlot(Tcell_seuratObj, group.by = 'cell_types_with_t_granular', cols = polychrome_pal)
DimPlot(Tcell_seuratObj, group.by = 'cell_types_with_t_coarse_v2', cols = polychrome_pal)

# dittoseq stacked barplot of T cell subset proportions
dittoBarPlot(
  object = Tcell_seuratObj,
  var = "cell_types_with_t_coarse_v2",
  group.by = "orig.ident")

dittoBarPlot(
  object = Tcell_seuratObj,
  var = "cell_types_with_t_granular",
  group.by = "orig.ident")

# CD8 T cells - add module score
CD8_seuratObj <- AddModuleScore(CD8_seuratObj,
                                features = list(CytotoxicityScore),
                                name="CytotoxicityScore")
VlnPlot(CD8_seuratObj, features = "CytotoxicityScore1", group.by = "orig.ident") + stat_summary(fun.y = median, geom='point', size = 15, colour = "white", shape = 95) 
