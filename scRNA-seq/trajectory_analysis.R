library(tidyverse)
library(Seurat)
library(patchwork)
library(scCustomize)
library(SeuratWrappers)
library(monocle3)
set.seed(1212)
options(scipen = 999)


# Figure 2G
# re-integrated epithelial populations ----

epithelial_seurat <- readRDS("../epithelial_integrated.RDS")
DimPlot_scCustom(epithelial_seurat, group.by = "integrated_snn_res.0.2")
DimPlot_scCustom(epithelial_seurat, group.by = "celltypes_granular")
exp_mat <- GetAssayData(epithelial_seurat, assay = "RNA", slot = "counts")
cell_meta <- epithelial_seurat[[]]
gene_anno <- as.data.frame(rownames(epithelial_seurat))
colnames(gene_anno) <- "gene_short_name"
rownames(gene_anno) <- rownames(epithelial_seurat)

cds <- new_cell_data_set(exp_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_anno)
cds <- make_cds(cds, epithelial_seurat)

cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds)
cds <- learn_graph(cds, use_partition = TRUE)
plot_cells(cds)

# pseudotime ----
panel_order <- c("Trp63", "Krt5", "Krt14", "Krt6a", "Krt8")

# choosing luminal cells as the starting point ----
cds_lum <- order_cells(cds)
colData(cds_lum)$pseudoime <- pseudotime(cds_lum)
plot_cells(cds_lum, color_cells_by = "pseudotime")
mini_lum <- cds_lum[rowData(cds_lum)$gene_short_name %in% 
              panel_order]
mini_lum_list <- make_minis_list(mini_lum)

cds_og_lum <- cds_og
colData(cds_og_lum)$pseudotime <- pseudotime(cds_lum)

pdf("trajectory_analysis_monocle3_integrated_true_granular_lum.pdf")
plot_cells(cds_lum, color_cells_by = "pseudotime", cell_size = 0.5)
plot_genes_in_pseudotime(mini_lum, color_cells_by = "celltypes_granular", panel_order = panel_order) +
  ggtitle("Joint pseudotime, colored by cell types") + 
  scale_color_manual(values = polychrome_pal)
plot_genes_in_pseudotime(mini_lum, color_cells_by = "orig.ident", panel_order = panel_order) + 
  ggtitle("Joint pseudotime, colored by samples") + 
  scale_color_manual(values = polychrome_pal)
plot_genes_in_pseudotime(mini_lum, color_cells_by = "integrated_snn_res.0.2", panel_order = panel_order) + 
  ggtitle("Joint pseudotime, colored by samples") + 
  scale_color_manual(values = polychrome_pal)
lapply(seq_along(mini_lum_list), plot_mini, mini = mini_lum_list, panel_order = panel_order, name = names(mini_lum_list))

dev.off()


# helper functions -------
make_minis_list <- function(mini) {
  mini_p18 <- mini[,colData(mini)$orig.ident == "P18"]
  mini_pf18 <- mini[,colData(mini)$orig.ident == "PF18"]
  mini_p12 <- mini[,colData(mini)$orig.ident == "P12"]
  mini_pf12 <- mini[,colData(mini)$orig.ident == "PF12"]
  mini_xp <- mini[,colData(mini)$orig.ident == "XP"]
  mini_xpf <- mini[,colData(mini)$orig.ident == "XPF"]
  list_mini <- list(mini_p18, mini_pf18, mini_p12, mini_pf12, mini_xp, mini_xpf)
  names(list_mini) <- c("P18", "PF18", "P12", "PF12", "XP", "XPF")
  return(list_mini)
}

plot_mini <- function(i, mini = mini, group = "celltypes_granular", panel_order = NULL, name) {
  print(name[[i]])
  plot_genes_in_pseudotime(mini[[i]], color_cells_by = group, panel_order = panel_order) +
    ggtitle(paste0("Joint pseudotime, plotting only ", name[[i]])) +
    scale_color_manual(values = polychrome_pal) + facet_wrap(~feature_label, scales = "fixed", ncol = 1, nrow =5)
}

make_cds <- function(cds, epithelial) {
  reducedDim(cds, type = "PCA") <- epithelial@reductions$pca@cell.embeddings
  cds@preprocess_aux$prop_var_expl <- epithelial@reductions$pca@stdev
  cds@int_colData@listData$reducedDims$UMAP <- epithelial@reductions$umap@cell.embeddings
  cds@reduce_dim_aux@listData[["UMAP"]] <-epithelial@reductions[["umap"]]@cell.embeddings
  cds@preprocess_aux$gene_loadings <- epithelial@reductions[["pca"]]@feature.loadings
  cds@clusters$UMAP_so$clusters <- epithelial@meta.data$your_clusters
  return(cds)
}
