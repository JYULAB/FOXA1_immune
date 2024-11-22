ibrary(tidyverse)
library(Seurat)
library(patchwork)
library(scCustomize)
set.seed(1212)
options(scipen = 999)

# each sample ----------
# P12
data_p12 <- Read10X("../../scRNA/full/M1385/filtered_feature_bc_matrix/")
p12 <- CreateSeuratObject(data_p12, project = "P12", min.cells = 3, min.features = 200)

p12[["percent.mt"]] <- PercentageFeatureSet(p12, pattern = "^mt-")
v <- VlnPlot(p12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0, combine = FALSE)
cutoff <- c(10000, 150000, 15)
v1 <- lapply(seq_along(v), function(x) v[[x]] + geom_hline(yintercept = cutoff) +
               theme(legend.position="none"))
v1[[1]] + v1[[2]] + v1[[3]]

p12 <- subset(p12, subset = nCount_RNA < 150000 & percent.mt < 15)
p12 <- NormalizeData(p12)
saveRDS(p12, "P12.RDS")

# PF12
data_pf12 <- Read10X("../../scRNA/full/M1386/filtered_feature_bc_matrix/")
pf12 <- CreateSeuratObject(data_pf12, project = "PF12", min.cells = 3, min.features = 200)

pf12[["percent.mt"]] <- PercentageFeatureSet(pf12, pattern = "^mt-")
v <- VlnPlot(pf12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0, combine = FALSE)
cutoff <- c(10000, 125000, 15)
v1 <- lapply(seq_along(v), function(x) v[[x]] + geom_hline(yintercept = cutoff) +
               theme(legend.position="none"))
v1[[1]] + v1[[2]] + v1[[3]]

pf12 <- subset(pf12, subset =  nCount_RNA < 125000 & percent.mt < 15)
pf12 <- NormalizeData(pf12)
saveRDS(pf12, "PF12.RDS")

# P18
data_p18 <- Read10X("../../scRNA/full/M1387/filtered_feature_bc_matrix/")
p18 <- CreateSeuratObject(data_p18, project = "P18", min.cells = 3, min.features = 200)
p18[["percent.mt"]] <- PercentageFeatureSet(p18, pattern = "^mt-")
v <- VlnPlot(p18, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0, combine = FALSE)
cutoff <- c(10000, 100000, 15)
v1 <- lapply(seq_along(v), function(x) v[[x]] + geom_hline(yintercept = cutoff) +
               theme(legend.position="none"))
v1[[1]] + v1[[2]] + v1[[3]]
p18 <- subset(p18, subset = nCount_RNA < 100000 & percent.mt < 15)
p18 <- NormalizeData(p18)
saveRDS(p18, "P18.RDS")

# PF18
data_pf18 <- Read10X("../../scRNA/full/M1418/filtered_feature_bc_matrix/")
pf18 <- CreateSeuratObject(data_pf18, project = "PF18", min.cells = 3, min.features = 200)
pf18[["percent.mt"]] <- PercentageFeatureSet(pf18, pattern = "^mt-")
v <- VlnPlot(pf18, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0, combine = FALSE)
cutoff <- c(10000, 100000, 15)
v1 <- lapply(seq_along(v), function(x) v[[x]] + geom_hline(yintercept = cutoff) +
               theme(legend.position="none"))
v1[[1]] + v1[[2]] + v1[[3]]

pf18 <- subset(pf18, subset =  nCount_RNA < 100000 & percent.mt < 15)
pf18 <- NormalizeData(pf18)
saveRDS(pf18, "PF18.RDS")

# XP
data_xp <- Read10X("../../scRNA/full/M1415/filtered_feature_bc_matrix/")
xp <- CreateSeuratObject(data_xp, project = "XP", min.cells = 3, min.features = 200)

xp[["percent.mt"]] <- PercentageFeatureSet(xp, pattern = "^mt-")
v <- VlnPlot(xp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0, combine = FALSE)
cutoff <- c(10000, 150000, 15)
v1 <- lapply(seq_along(v), function(x) v[[x]] + geom_hline(yintercept = cutoff) +
               theme(legend.position="none"))
v1[[1]] + v1[[2]] + v1[[3]]

xp <- subset(xp, subset = nCount_RNA < 100000 & percent.mt < 15)
xp <- NormalizeData(xp)
saveRDS(xp, "XP.RDS")

# XPF

data_xpf <- Read10X("../../scRNA/full/M1475/filtered_feature_bc_matrix/")
xpf <- CreateSeuratObject(data_xpf, project = "XPF", min.cells = 3, min.features = 200)

xpf[["percent.mt"]] <- PercentageFeatureSet(xpf, pattern = "^mt-")
v <- VlnPlot(xpf, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0, combine = FALSE)
v1 <- lapply(seq_along(v), function(x) v[[x]] + geom_hline(yintercept = cutoff) +
               theme(legend.position="none"))
v1[[1]] + v1[[2]] + v1[[3]]

xpf <- subset(xpf, subset = nCount_RNA < 100000 & percent.mt < 15)
xpf <- NormalizeData(xpf)
saveRDS(xpf, "XPF.RDS")