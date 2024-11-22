library(Seurat)
library(tidyverse)
library(patchwork)
library(scCustomize)
library(gridExtra)
set.seed(1212)
options(scipen = 999)

polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]

# processing each sample -------
data_pf15 <- Read10X("../../pip-seq/PF15/sensitivity_5/")
pf15 <- CreateSeuratObject(data_pf15, project = "PF15", min.cells = 3, min.features = 200)

pf15[["percent.mt"]] <- PercentageFeatureSet(pf15, pattern = "^mt-")
v <- VlnPlot(pf15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0, combine = FALSE)
cutoff <- c(4000, 10000, 10)
v1 <- lapply(seq_along(v), function(x) v[[x]] + geom_hline(yintercept = cutoff[x]) +
               theme(legend.position="none"))
v_tot <- v1[[1]] + v1[[2]] + v1[[3]] + plot_annotation("PF15")
v_tot
pf15 <- subset(pf15, subset = nCount_RNA < 10000 & nFeature_RNA < 4000 & percent.mt < 10)
pf15 <- NormalizeData(pf15)

data_p15 <- Read10X("../../pip-seq/M1488/sensitivity_5/")
p15 <- CreateSeuratObject(data_p15, project = "P15", min.cells = 3, min.features = 200)

p15[["percent.mt"]] <- PercentageFeatureSet(p15, pattern = "^mt-")
w <- VlnPlot(p15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0, combine = FALSE)
cutoff <- c(4000, 10000, 10)
w1 <- lapply(seq_along(w), function(x) w[[x]] + geom_hline(yintercept = cutoff[x]) +
               theme(legend.position="none"))
w_tot <- w1[[1]] + w1[[2]] + w1[[3]]  + plot_annotation("P15")
w_tot
v_tot / w_tot

p15 <- subset(p15, subset = nCount_RNA < 10000 & nFeature_RNA < 4000 & percent.mt < 10)
p15 <- NormalizeData(p15)

# integrating ------
samples_list <- list(p15, pf15)
samples_list <- lapply(X = samples_list, FUN = SCTransform, vst.flavor = "v2")
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)
samples_list <- PrepSCTIntegration(object.list = samples_list, anchor.features = features)

anchors_cca <- FindIntegrationAnchors(object.list = samples_list, normalization.method = "SCT",
                                      anchor.features = features)
integrated <- IntegrateData(anchorset = anchors_cca, normalization.method = "SCT")
DefaultAssay(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)

DefaultAssay(integrated) <- "integrated"
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.4)
DimPlot_scCustom(integrated, colors_use = polychrome_pal)

# cell typing with markers ----
Idents(integrated) <- "integrated_snn_res.0.4"
cell_types <- c("0" = "Macrophage", "1" = "B Cell", "2" = "T Cell (NK/CD8+)", "3" = "Th17", "4" = "M2 Macrophage", 
                "5" = "T Cell (CD4+)", "6" = "Dendritic Cell", "7" = "T Cell (NK/CD8+)", "8"= "T Cell (Foxp3+)", 
                "9" = "Tumor Epithelial", "10" = "Cell Cycle CD8 T Cell", "11" = "Monocyte", "12" = "T Cell", 
                "13" = "B Cell", "14" = "Macrophage")
Idents(integrated) <- "integrated_snn_res.0.4"
integrated <- RenameIdents(integrated, cell_types)
integrated$cell_types <- Idents(integrated)
DimPlot_scCustom(integrated, group.by = "cell_types", colors_use = polychrome_pal)
integrated$cell_types <- factor(integrated$cell_types, 
                                levels = c("B Cell", "Dendritic Cell", "Macrophage", "M2 Macrophage", "Monocyte",
                                           "T Cell", "T Cell (CD4+)", "T Cell (NK/CD8+)", "T Cell (Foxp3+)", "Th17",
                                           "Cell Cycle CD8 T Cell","Tumor Epithelial"))
Idents(integrated) <- "cell_types"

DefaultAssay(integrated) <- "RNA"
pdf("reseq_integrated_pip-seq_gene_expression.pdf")
DimPlot_scCustom(integrated, colors_use = polychrome_pal, group.by = "orig.ident") + ggtitle("Integrated")
VlnPlot_scCustom(integrated, c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
FeaturePlot_scCustom(integrated, c("nFeature_RNA", "nCount_RNA", "percent.mt"))

FeaturePlot_scCustom(integrated, c("Ptprc", "Itgam", "Ly6g", "Ms4a1"))
FeaturePlot_scCustom(integrated, c("S100a8", "S100a9", "Cd14")) + plot_annotation(title = "Monocytes") +
  plot_layout(ncol = 2)
FeaturePlot_scCustom(integrated, c("Adgre1", "Mafb", "Lyz2", "C1qb")) + plot_annotation(title = "Macrophages")
FeaturePlot_scCustom(integrated, c("Mrc1", "Cd163", "Arg1")) + plot_annotation(title = "M2 Macrophages") +
  plot_layout(ncol = 2)
FeaturePlot_scCustom(integrated, c("Cd86", "Cd80", "Nos2")) + plot_annotation(title = "M1 Macrophages") +
  plot_layout(ncol = 2)

FeaturePlot_scCustom(integrated, c("Cd3e", "Cd4", "Cd8a", "Foxp3")) + plot_annotation(title = "T cells")
FeaturePlot_scCustom(integrated, c("Havcr2", "Tigit")) + plot_annotation(title = "Exhausted T Cells") + plot_layout(nrow = 2)
FeaturePlot_scCustom(integrated, c("Gzmb", "Prf1")) + plot_annotation(title = "Cytotoxic T cells") + plot_layout(nrow = 2)

FeaturePlot_scCustom(integrated, c("Klrd1", "Nkg7")) + plot_annotation(title = "NK cells") + 
  plot_layout(ncol= 2, nrow =2)
FeaturePlot_scCustom(integrated, c("Itgax", "H2-DMb2", "Xcr1", "Batf3"))  + 
  plot_annotation(title = "DCs")
FeaturePlot_scCustom(integrated, c("Cd209a", "Flt3"))  + 
  plot_annotation(title = "DCs") +
  plot_layout(nrow = 2)

FeaturePlot_scCustom(integrated, c("Pdcd1", "Cd274", "Ctla4", "Lag3")) + plot_annotation(title = "Immune Checkpoints")
FeaturePlot_scCustom(integrated, c("Epcam", "Mki67")) + plot_layout(nrow = 2)
dev.off()


# subsetting on T cells and annotating in more details -------
Idents(integrated) <- "integrated_snn_res.0.4"
DimPlot_scCustom(integrated)
t_cell_og <- subset(integrated, idents = c("2", "3", "5", "7", "8", "10", "12"))
DimPlot_scCustom(t_cell_og)
t_cell_list <- SplitObject(t_cell_og, split.by = "orig.ident")
t_cell_list <- lapply(X = t_cell_list, FUN = SCTransform, vst.flavor = "v2")
t_cell_features <- SelectIntegrationFeatures(object.list = t_cell_list, nfeatures = 3000)
t_cell_list <- PrepSCTIntegration(object.list = t_cell_list, anchor.features = t_cell_features)

t_cell_anchors_cca <- FindIntegrationAnchors(object.list = t_cell_list, normalization.method = "SCT",
                                             anchor.features = t_cell_features)
t_cell_integrated <- IntegrateData(anchorset = t_cell_anchors_cca, normalization.method = "SCT")
t_cell_integrated <- RunPCA(t_cell_integrated)
t_cell_integrated <- RunUMAP(t_cell_integrated, dims = 1:30)
t_cell_integrated$integrated_snn_res.0.4_og <- t_cell_integrated$integrated_snn_res.0.4

DefaultAssay(t_cell_integrated) <- "integrated"
t_cell_integrated <- FindNeighbors(t_cell_integrated, dims = 1:30)
t_cell_integrated <- FindClusters(t_cell_integrated, resolution = 0.6)

# making cluster 5 separated based on CD3e expression ------
Idents(t_cell_integrated) <- "integrated_snn_res.0.6"
cluster5_cd3e_neg <- WhichCells(t_cell_integrated, expression = Cd3e == 0, idents = "5")
Cell_Highlight_Plot(t_cell_integrated, cells_highlight = list(Cd3e_neg = cluster5_cd3e_neg), highlight_color = c("forestgreen"))
t_cell_integrated$manual_clusters <- t_cell_integrated$integrated_snn_res.0.6
levels(t_cell_integrated$manual_clusters) <- c(levels(t_cell_integrated$manual_clusters), "13") # assigning new cluster 13
t_cell_integrated$manual_clusters[names(t_cell_integrated$manual_clusters) %in% cluster5_cd3e_neg] <- "13"

cluster9_nkg7_neg <- WhichCells(t_cell_integrated, expression = Nkg7 == 0, idents = "9")
Cell_Highlight_Plot(t_cell_integrated, cells_highlight = list(Nkg7_neg = cluster9_nkg7_neg), highlight_color = c("forestgreen"))
levels(t_cell_integrated$manual_clusters) <- c(levels(t_cell_integrated$manual_clusters), "14") # assigning new cluster 14
t_cell_integrated$manual_clusters[names(t_cell_integrated$manual_clusters) %in% cluster9_nkg7_neg] <- "14" 

# t cell type granular ----
t_cell_types_granular <- c("0" = "Naive/Memory T helper", "1" = "Th17", "2" = "CTL1", 
                           "3" = "CTL2", "4" = "Th17", 
                           "5" = "CTL3", "6" = "Tregs", 
                           "7" = "Cycling CD8","8"= "Tregs",
                           "9" = "NK", "10" = "Naive/Memory T helper",
                           "11" = "CTL1", "12" = "CTL1", "13" = "NK", "14" = "CTL2")
Idents(t_cell_integrated) <- "manual_clusters"
t_cell_integrated <- RenameIdents(t_cell_integrated, t_cell_types_granular)
t_cell_integrated$t_cell_types_granular <- Idents(t_cell_integrated)
DimPlot_scCustom(t_cell_integrated, group.by = "t_cell_types_granular", colors_use = polychrome_pal, pt.size = 1)
Idents(t_cell_integrated) <- "t_cell_types_granular"


# t cell type coarse ----
t_cell_types_coarse <- c("0" = "T helper", "1" = "T helper", "2" = "CD8 T Cell", 
                         "3" = "CD8 T Cell", "4" = "T helper", 
                         "5" = "CD8 T Cell", "6" = "T helper", 
                         "7" = "CD8 T Cell", "8"= "T helper",
                         "9" = "NK", "10" = "T helper",
                         "11" = "CD8 T Cell", "12" = "CD8 T Cell", "13" = "NK", "14" = "CD8 T Cell")
Idents(t_cell_integrated) <- "manual_clusters"
t_cell_integrated <- RenameIdents(t_cell_integrated, t_cell_types_coarse)
t_cell_integrated$t_cell_types_coarse <- Idents(t_cell_integrated)

pdf("reseq_t_cells.pdf")
DimPlot_scCustom(t_cell_og, colors_use = polychrome_pal) + plot_annotation("T Cell Subset")
DimPlot_scCustom(t_cell_integrated, colors_use = polychrome_pal) + plot_annotation("T Cell Subset")
DefaultAssay(t_cell_integrated) <- "RNA"
FeaturePlot_scCustom(t_cell_integrated, c("Cd3e", "Cd4", "Cd8a", "Foxp3")) + plot_annotation(title = "T cells")
FeaturePlot_scCustom(t_cell_integrated, c("Havcr2", "Tigit")) + plot_annotation(title = "Exhausted T Cells") + plot_layout(nrow = 2)
FeaturePlot_scCustom(t_cell_integrated, c("Gzmb", "Prf1")) + plot_annotation(title = "Cytotoxic T cells") + plot_layout(nrow = 2)
FeaturePlot_scCustom(t_cell_integrated, c("Klrd1", "Nkg7")) + plot_annotation(title = "NK cells") + 
  plot_layout(ncol= 2, nrow =2)

FeaturePlot_scCustom(t_cell_integrated, c("Gzma", "Top2a", "Pdcd1", "Cdk1"))
FeaturePlot_scCustom(t_cell_integrated, c("Trbc2", "Cd28", "Icos", "Thy1"))
FeaturePlot_scCustom(t_cell_integrated, c("Il7r", "Fcgr4", "Ppbp", "Cst3"))

dev.off()

# adding back to entire object ------
t_cell_in_main <- Cells(integrated) %in% names(t_cell_integrated$cell_types)

integrated$cell_types_with_t_coarse <- integrated$cell_types
integrated$cell_types_with_t_coarse <- as.character(integrated$cell_types_with_t_coarse)
DimPlot_scCustom(integrated, group.by = "cell_types_with_t_coarse", colors_use = polychrome_pal)
integrated$cell_types_with_t_coarse[t_cell_in_main] <- as.character(t_cell_integrated$t_cell_types_coarse)
DimPlot_scCustom(integrated, group.by = "cell_types_with_t_coarse", colors_use = polychrome_pal)

integrated$cell_types_with_t_granular <- integrated$cell_types
integrated$cell_types_with_t_granular <- as.character(integrated$cell_types_with_t_granular)
DimPlot_scCustom(integrated, group.by = "cell_types_with_t_granular", colors_use = polychrome_pal)
integrated$cell_types_with_t_granular[t_cell_in_main] <- as.character(t_cell_integrated$t_cell_types_granular)
DimPlot_scCustom(integrated, group.by = "cell_types_with_t_granular", colors_use = polychrome_pal)
DefaultAssay(integrated) <- "RNA"
saveRDS(integrated, "pip-seq_integrated.RDS")


