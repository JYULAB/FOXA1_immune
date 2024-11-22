library(tidyverse)
library(Seurat)
library(patchwork)
library(scCustomize)
set.seed(1212)
options(scipen = 999)
polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]
stepped <- DiscretePalette_scCustomize(num_colors = 24, palette = "stepped")


p12 <- readRDS("../P12.RDS")
pf12 <- readRDS("../PF12.RDS")
p18 <- readRDS("../P18.RDS")
pf18 <- readRDS("../PF18.RDS")
xp <- readRDS("../XP.RDS")
xpf <- readRDS("../XPF.RDS")

# merge samples -------

samples_list <- list(p12, pf12, p18, pf18, xp, xpf)
samples_list <- lapply(X = samples_list, FUN = SCTransform, vst.flavor = "v2")

# integrations -----

features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)
samples_list <- PrepSCTIntegration(object.list = samples_list, anchor.features = features)

anchors_cca <- FindIntegrationAnchors(object.list = samples_list, normalization.method = "SCT",
                                      anchor.features = features)
integrated <- IntegrateData(anchorset = anchors_cca, normalization.method = "SCT")
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)
DimPlot_scCustom(integrated, group.by = "orig.ident", colors_use = polychrome_pal) + ggtitle("Integrated")

integrated$timepoint <- integrated$orig.ident %>% str_extract("\\d{2}")
integrated$timepoint <- ifelse(is.na(integrated$timepoint), "X", integrated$timepoint)
DimPlot_scCustom(integrated, colors_use = polychrome_pal, group.by = "timepoint") + 
  labs(title = "Timepoint")

integrated$genotype <- integrated$orig.ident %>% str_extract("[:alpha:]*")
DimPlot_scCustom(integrated, colors_use = polychrome_pal, group.by = "genotype") + 
  labs(title = "Genotype")
saveRDS(integrated, "integrated.RDS")

# plotting cell type markers -----

pdf("integrated_gene_expression.pdf")
DefaultAssay(integrated) <- "RNA"
FeaturePlot_scCustom(integrated, c("Ar", "Krt8", "Cd24a", "Krt18")) + plot_annotation(title = "Luminal")
FeaturePlot_scCustom(integrated, c("Spink1", "Pbsn")) + plot_annotation(title = "Additional Luminal") +
  plot_layout(nrow = 2)
FeaturePlot_scCustom(integrated, c("Trp63", "Krt5", "Krt14", "Lgals7")) + plot_annotation(title = "Basal")
FeaturePlot_scCustom(integrated, c("Pate4", "Pax2", "Svs2")) + plot_annotation(title = "Seminal vesicles")
FeaturePlot_scCustom(integrated, c("Ascl1", "Chga", "Foxa2", "Eno2")) + plot_annotation(title = "NE")
FeaturePlot_scCustom(integrated, c("Dsc2", "Dsg2", "Krt13", "Krt16")) + plot_annotation(title = "Squamous")
FeaturePlot_scCustom(integrated, c("Dcn", "Lum", "Vim")) + plot_annotation(title = "Fibroblast")
FeaturePlot_scCustom(integrated, c("Agr2", "Krt7", "Tff3")) + plot_annotation(title = "Club Epithelial") + 
  plot_layout(nrow = 2)
FeaturePlot_scCustom(integrated, c("Krt13", "Krt19", "Cldn4")) + plot_annotation(title = "Hillock Epithelial")
FeaturePlot_scCustom(integrated, c("Ackr1", "Cldn5")) + plot_annotation(title = "Endothelial") + 
  plot_layout(nrow = 2)

FeaturePlot_scCustom(integrated, c("Cd14", "S100a8", "S100a9", "Cd68")) + plot_annotation(title = "Monocytes")
FeaturePlot_scCustom(integrated, c("Ptprc", "Itgam", "Adgre1", "Itgax")) + plot_annotation(title = "Myeloid cells")
FeaturePlot_scCustom(integrated, c("Cd3e", "Ms4a1", "Klrb1c")) + plot_annotation(title = "Lymphoid Cells")
dev.off()

# annotation ------

DefaultAssay(integrated) <- "integrated"
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.2)

celltypes_main <- c("0" = "Luminal", "1" = "Basal", "2" = "Basal", "3" = "Immune",
                      "4" = "Luminal", "5" = "Basal", "6" = "Seminal Vesicles", "7" = "Immune",
                      "8" = "Luminal", "9" = "Basal",
                      "10" = "Stroma", "11" = "Luminal", "12" = "Immune", "13" = "Seminal Vesicles",
                      "14" = "Immune", "15" = "Basal", "16" = "Stroma")
Idents(integrated) <- "integrated_snn_res.0.2"
integrated <- RenameIdents(integrated, celltypes_main)
integrated$celltypes_main <- factor(Idents(integrated), 
                                    levels = c("Luminal", "Basal", "Seminal Vesicles",
                                               "Stroma", "Immune"))
Idents(integrated) <- "celltypes_main"
DimPlot_scCustom(integrated, colors_use = polychrome_pal)

celltypes_coarse <- c("0" = "Luminal", "1" = "Basal", "2" = "Basal", "3" = "Monocytes",
                      "4" = "Luminal", "5" = "Basal", "6" = "Seminal Vesicles", "7" = "Macrophage",
                      "8" = "Luminal", "9" = "Basal",
                      "10" = "Stroma", "11" = "Luminal", "12" = "T Cells", "13" = "Seminal Vesicles",
                      "14" = "B Cells", "15" = "Basal", "16" = "Stroma")
Idents(integrated) <- "integrated_snn_res.0.2"
integrated <- RenameIdents(integrated, celltypes_coarse)
integrated$celltypes_coarse <- factor(Idents(integrated), 
                                      levels = c("Luminal", "Basal", "Seminal Vesicles",
                                                 "Stroma", "Spink1", "Monocytes", "Macrophage",
                                                 "B Cells", "T Cells"))
Idents(integrated) <- "celltypes_coarse"
DimPlot_scCustom(integrated, colors_use = polychrome_pal)

celltypes_granular <- c("0" = "Luminal", "1" = "Basal", "2" = "Basal Squamous", "3" = "Monocytes",
                      "4" = "Spink1", "5" = "Basal", "6" = "Seminal Vesicles", "7" = "Macrophage",
                      "8" = "Pbsn", "9" = "Basal Proliferating",
                      "10" = "Fibroblasts", "11" = "Club Epithelial", "12" = "T cells", "13" = "Seminal Vesicles",
                      "14" = "B cells", "15" = "Basal", "16" ="Endothelial")
Idents(integrated) <- "integrated_snn_res.0.2"
integrated <- RenameIdents(integrated, celltypes_granular)
integrated$celltypes_granular <- Idents(integrated)
integrated$celltypes_granular <- factor(integrated$celltypes_granular, levels = 
                                          c("Luminal", "Spink1", "Pbsn", 
                                           "Club Epithelial", "Seminal Vesicles",
                                           "Basal", "Basal Squamous", "Basal Proliferating",
                                           "Fibroblasts", "Endothelial",
                                           "Monocytes", "Macrophage", "B cells", "T cells"))

Idents(integrated) <- "celltypes_granular"
DimPlot_scCustom(integrated, colors_use = polychrome_pal)


# selecting for epithelial cells and re-integrating ------

library(monocle3)
library(SeuratWrappers)

DimPlot_scCustom(integrated, colors_use = polychrome_pal, group.by = "celltypes_coarse")
Idents(integrated) <- "celltypes_main" 
epithelial <- subset(integrated, idents = c("Luminal", "Basal", "Spink1"))
DimPlot_scCustom(epithelial, group.by = "integrated_snn_res.0.2")
DimPlot_scCustom(epithelial, group.by = "celltypes_main")
Meta_Highlight_Plot(epithelial, meta_data_column = "integrated_snn_res.0.2", 
                    meta_data_highlight = "5")
DefaultAssay(epithelial) <- "RNA"
saveRDS(epithelial, "epithelial_correct.RDS")

# re-integrating epithelial subsets -----

#epithelial <- readRDS("../epithelial_correct.RDS")
DimPlot_scCustom(epithelial, group.by = "celltypes_granular")
epithelial_list <- SplitObject(epithelial, split.by = "orig.ident")
epithelial_list <- lapply(X = epithelial_list, FUN = SCTransform, vst.flavor = "v2")
epithelial_features <- SelectIntegrationFeatures(object.list = epithelial_list, nfeatures = 3000)
epithelial_list <- PrepSCTIntegration(object.list = epithelial_list, anchor.features = epithelial_features)

epithelial_anchors_cca <- FindIntegrationAnchors(object.list = epithelial_list, normalization.method = "SCT",
                                                 anchor.features = epithelial_features)
epithelial_integrated <- IntegrateData(anchorset = epithelial_anchors_cca, normalization.method = "SCT")
epithelial_integrated <- RunPCA(epithelial_integrated)
epithelial_integrated <- RunUMAP(epithelial_integrated, dims = 1:30)

epithelial_integrated$integrated_snn_res.0.2_og <- epithelial_integrated$integrated_snn_res.0.2
epithelial_integrated <- FindNeighbors(epithelial_integrated, dims = 1:30)
epithelial_integrated <- FindClusters(epithelial_integrated, resolution = 0.2)
DimPlot_scCustom(epithelial_integrated, group.by = "integrated_snn_res.0.2")
DefaultAssay(epithelial_integrated) <- "RNA"
saveRDS(epithelial_integrated, "epithelial_integrated.RDS")

# slim down seurat object for sharing ----
diet_seurat <- DietSeurat(integrated, assays = "RNA", scale.data = TRUE, 
                          dimreducs = c("pca", "umap"))
saveRDS(diet_seurat, "minimal_seurat.RDS")
