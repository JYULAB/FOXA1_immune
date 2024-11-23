# Analysis of previously published human PCa scRNA-seq data - Hirz et al. Nat Comms.  2023
# Citation - data source - GSE181294
# Hirz, T., Mei, S., Sarkar, H. et al. Dissecting the immune suppressive human prostate tumor microenvironment 
# via integrated single-cell and spatial transcriptomic analyses. Nat Commun 14, 663 (2023). 
# https://doi.org/10.1038/s41467-023-36325-2

# Fig 7

# Load Packages:
library(Seurat)
library(SeuratData)
library(scCustomize)
library(ggplot2)
library(dittoSeq)
library(dplyr) 
library(viridis)

# Read in seurat object
seuratObj = readRDS("hirz_et_al_seurat_object.rds")
seuratObj = UpdateSeuratObject(object = seuratObj)

# color palettes
varibow_pal <- scCustomize::DiscretePalette_scCustomize(num_colors = 50, palette = "varibow")
polychrome_pal <- scCustomize::DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
polychrome_pal <- polychrome_pal[c(1, 3, 6:36)]

# dimplot
DimPlot(seuratObj, group.by = "orig.ident", cols = varibow_pal, shuffle = TRUE)
DimPlot(seuratObj, group.by = "status", cols = polychrome_pal, shuffle = TRUE)

# Replicate Hirz 2023 cell type UMAP
Epithelial <- c("Epitheial","Tumor")
Thelper <- c("CD4")
Treg <- c("Treg")
CTL <- c("CTL")
CyclingT <- c("Cycling T")
Bcell <- c("B cells")
Plasma <- c("Plasma cells")
NK <- c("NK")
Macrophage <- c("Macrophage")
Monocytes <- c("Monocytes")
Mast <- c("Mast cells")
mDC <- c("mDC")
pDC <- c("PDC")
Pericytes <- c("Pericytes")
Fibroblasts <- c("Fibroblasts")
Endothelial <- c("Endothelial cell")

cl.mat <- cbind(c(rep("Epithelial", length(Epithelial)), rep("T helper", length(Thelper)), 
                  rep("Treg", length(Treg)), rep("CTL", length(CTL)), 
                  rep("Cycling T", length(CyclingT)), rep("B cell", length(Bcell)), 
                  rep("Plasma cell", length(Plasma)), rep("NK", length(NK)), 
                  rep("Macrophage", length(Macrophage)), rep("Monocytes", length(Monocytes)),
                  rep("Mast cell", length(Mast)), rep("mDC", length(mDC)), rep("pDC", length(pDC)), 
                  rep("Pericytes", length(Pericytes)), rep("Fibroblasts", length(Fibroblasts)), 
                  rep("Endothelial", length(Endothelial))
),
c(Epithelial, Thelper, Treg, CTL, CyclingT, Bcell, Plasma, NK, Macrophage, Monocytes, Mast, mDC, pDC, Pericytes, Fibroblasts,  Endothelial))

cl.vec <- seuratObj$cells1
ct.vec <- rep(NA, length(cl.vec))
for(x in unique(cl.mat[,1])){
  cl.x <- cl.mat[cl.mat[,1]==x,2]
  ct.vec[which(cl.vec%in%cl.x)] <- x
}
seuratObj$cells3 <- ct.vec
seuratObj$cells3 <- factor(seuratObj$cells3, 
                                     levels = c('Epithelial', 'T helper','Treg','CTL',
                                                'Cycling T','B cell','Plasma cell',
                                                'NK','Macrophage','Monocytes','Mast cell',
                                                'mDC','pDC','Pericytes','Fibroblasts','Endothelial'))

DimPlot(seuratObj, group.by = 'cells3', cols = polychrome_pal, shuffle = TRUE)

#-------------------------------------------------------------------------------
# subset on tumor samples
Tumor_seuratObj <- seuratObj[,seuratObj$status == "Tumor"]
DimPlot(Tumor_seuratObj, group.by = "orig.ident", cols = polychrome_pal, shuffle = TRUE)

# remove sample SCG-PCA5-T-LG because too few epithelial cells
Tumor_seuratObj <- Tumor_seuratObj[,Tumor_seuratObj$orig.ident != "SCG-PCA5-T-LG"]
DimPlot(Tumor_seuratObj, group.by = "orig.ident", cols = polychrome_pal, shuffle = TRUE)

# order samples from FOXA1 high to low
Idents(Tumor_seuratObj) <- Tumor_seuratObj$cells3
VlnPlot(Tumor_seuratObj, features = "FOXA1", idents = "Epithelial", group.by = "orig.ident", sort = 'increasing')+stat_summary(fun.y = mean, geom='point', size = 10, colour = "red", shape = 95)
Tumor_seuratObj$orig.ident <- factor(Tumor_seuratObj$orig.ident, levels = c('SCG-PCA17-T-LG',	
                                                                            'SCG-PCA9-T-LG',	'SCG-PCA12-T-LG',	'SCG-PCA18-T-LG',	'SCG-PCA3-T-LG',	
                                                                            'SCG-PCA11-T-LG',	'SCG-PCA6-T-HG',	'SCG-PCA15-T-HG',	'SCG-PCA16-T-HG',	
                                                                            'SCG-PCA7-T-HG',	'SCG-PCA8-T-HG',	'SCG-PCA19-T-HG',	'SCG-PCA21-T-LG',	
                                                                            'SCG-PCA4-T-HG',	'SCG-PCA10-T-LG',	'SCG-PCA20-T-LG',	'SCG-PCA22-T-HG',	
                                                                            'SCG-PCA5-T-LG'))

# Add FOXA1 low vs high groups for tumor samples
FOXA1high <- c('SCG-PCA17-T-LG', 'SCG-PCA9-T-LG',	'SCG-PCA12-T-LG',	'SCG-PCA18-T-LG',	'SCG-PCA3-T-LG',	
               'SCG-PCA11-T-LG',	'SCG-PCA6-T-HG',	'SCG-PCA15-T-HG')
FOXA1low <- c('SCG-PCA16-T-HG', 'SCG-PCA7-T-HG',	'SCG-PCA8-T-HG',	'SCG-PCA19-T-HG',	'SCG-PCA21-T-LG',	
              'SCG-PCA4-T-HG',	'SCG-PCA10-T-LG',	'SCG-PCA20-T-LG',	'SCG-PCA22-T-HG',	'SCG-PCA5-T-LG')

cl.mat <- cbind(c(rep("FOXA1 high", length(FOXA1high)), rep("FOXA1 low", length(FOXA1low))
),
c(FOXA1high, FOXA1low))

cl.vec <- Tumor_seuratObj$orig.ident
ct.vec <- rep(NA, length(cl.vec))
for(x in unique(cl.mat[,1])){
  cl.x <- cl.mat[cl.mat[,1]==x,2]
  ct.vec[which(cl.vec%in%cl.x)] <- x
}
Tumor_seuratObj$FOXA1status <- ct.vec
Tumor_seuratObj$FOXA1status <- factor(Tumor_seuratObj$FOXA1status, 
                                levels = c('FOXA1 high', 'FOXA1 low'))
DimPlot(Tumor_seuratObj, group.by = 'orig.ident', split.by = 'FOXA1status', cols = varibow_pal, shuffle = TRUE)

#-------------------------------------------------------------------------------
# subset on epithelial cells in tumor samples
Tumor_Epi_seuratObj <- Tumor_seuratObj[,Tumor_seuratObj$cells3 == "Epithelial"]
DimPlot(Tumor_Epi_seuratObj, group.by = "orig.ident", cols = polychrome_pal, shuffle = TRUE)

# subset on Macrophages in tumor samples 
Mac_seuratObj <- Tumor_seuratObj[,Tumor_seuratObj$cells3 == "Macrophage"]
DimPlot(Mac_seuratObj, group.by = "orig.ident", cols = polychrome_pal, shuffle = TRUE)

# subset on CTL in tumor samples 
CTL_seuratObj <- Tumor_seuratObj[,Tumor_seuratObj$cells3 == "CTL"]
DimPlot(CTL_seuratObj, group.by = "orig.ident", cols = polychrome_pal, shuffle = TRUE)

#-------------------------------------------------------------------------------
# Violin plots

# Violin plot of FOXA1 expression in epitheial cells - Fig S7
VlnPlot(Tumor_Epi_seuratObj, features = "FOXA1",  group.by = "orig.ident")+stat_summary(fun.y = mean, geom='point', size = 10, colour = "red", shape = 95)

# Violin plots of lineage marker gene expression - Fig 7
VlnPlot(Tumor_Epi_seuratObj, features = "FOXA1", group.by = "FOXA1status")
VlnPlot(Tumor_Epi_seuratObj, features = c("KRT8","AR","KLK3","TP63","KRT5","KRT13"), group.by = "FOXA1status") 

#-------------------------------------------------------------------------------
# Signature gene lists
# MSigDB gene sets
TGFB_signaling_H <- c(	'ACVR1',	'APC',	'ARID4B',	'BCAR3',	'BMP2',	'BMPR1A',	'BMPR2',	'CDH1',	'CDK9',	'CDKN1C',	'CTNNB1',	
                       'ENG',	'FKBP1A',	'FNTA',	'FURIN',	'HDAC1',	'HIPK2',	'ID1',	'ID2',	'ID3',	'IFNGR2',	'JUNB',	'KLF10',	
                       'LEFTY2',	'LTBP2',	'MAP3K7',	'NCOR2',	'NOG',	'PMEPA1',	'PPM1A',	'PPP1CA',	'PPP1R15A',	'RAB31',	'RHOA',	
                       'SERPINE1',	'SKI',	'SKIL',	'SLC20A1',	'SMAD1',	'SMAD3',	'SMAD6',	'SMAD7',	'SMURF1',	'SMURF2',	'SPTBN1',	
                       'TGFB1',	'TGFBR1',	'TGIF1',	'THBS1',	'TJP1',	'TRIM33',	'UBE2D3',	'WWTR1',	'XIAP')																																																																																																																																																		
TNFA_signaling_H <- c(	'ABCA1',	'AREG',	'ATF3',	'ATP2B1',	'B4GALT1',	'B4GALT5',	'BCL2A1',	'BCL3',	'BCL6',	'BHLHE40',	'BIRC2',	'BIRC3',	
                       'BMP2',	'BTG1',	'BTG2',	'BTG3',	'CCL2',	'CCL20',	'CCL4',	'CCL5',	'CCND1',	'CCNL1',	'CCRL2',	'CD44',	'CD69',	'CD80',	'CD83',	
                       'CDKN1A',	'CEBPB',	'CEBPD',	'CFLAR',	'CLCF1',	'CSF1',	'CSF2',	'CXCL1',	'CXCL10',	'CXCL11',	'CXCL2',	'CXCL3',	'CXCL6',	
                       'CXCR7',	'CYR61',	'DDX58',	'DENND5A',	'DNAJB4',	'DRAM1',	'DUSP1',	'DUSP2',	'DUSP4',	'DUSP5',	'EDN1',	'EFNA1',	'EGR1',	
                       'EGR2',	'EGR3',	'EHD1',	'EIF1',	'ETS2',	'F2RL1',	'F3',	'FJX1',	'FOS',	'FOSB',	'FOSL1',	'FOSL2',	'FUT4',	'G0S2',	'GADD45A',	
                       'GADD45B',	'GCH1',	'GEM',	'GFPT2',	'GPR183',	'HBEGF',	'HES1',	'ICAM1',	'ICOSLG',	'ID2',	'IER2',	'IER3',	'IER5',	'IFIH1',	
                       'IFIT2',	'IFNGR2',	'IL12B',	'IL15RA',	'IL18',	'IL1A',	'IL1B',	'IL23A',	'IL6',	'IL6ST',	'IL7R',	'INHBA',	'IRF1',	'IRS2',	
                       'JAG1',	'JUN',	'JUNB',	'KDM6B',	'KLF10',	'KLF2',	'KLF4',	'KLF6',	'KLF9',	'KYNU',	'LAMB3',	'LDLR',	'LIF',	'LITAF',	'MAFF',	
                       'MAP2K3',	'MAP3K8',	'MARCKS',	'MCL1',	'MSC',	'MXD1',	'MYC',	'NAMPT',	'NFAT5',	'NFE2L2',	'NFIL3',	'NFKB1',	'NFKB2',	
                       'NFKBIA',	'NFKBIE',	'NINJ1',	'NR4A1',	'NR4A2',	'NR4A3',	'OLR1',	'PANX1',	'PDE4B',	'PDLIM5',	'PER1',	'PFKFB3',	'PHLDA1',	
                       'PHLDA2',	'PLAU',	'PLAUR',	'PLEK',	'PLK2',	'PMEPA1',	'PNRC1',	'PPAP2B',	'PPP1R15A',	'PTGER4',	'PTGS2',	'PTPRE',	'PTX3',	
                       'RCAN1',	'REL',	'RELA',	'RELB',	'RHOB',	'RIPK2',	'RNF19B',	'SAT1',	'SDC4',	'SERPINB2',	'SERPINB8',	'SERPINE1',	'SGK1',	'SIK1',	
                       'SLC16A6',	'SLC2A3',	'SLC2A6',	'SMAD3',	'SNN',	'SOCS3',	'SOD2',	'SPHK1',	'SPSB1',	'SQSTM1',	'STAT5A',	'TANK',	'TAP1',	'TGIF1',	
                       'TIPARP',	'TLR2',	'TNC',	'TNF',	'TNFAIP2',	'TNFAIP3',	'TNFAIP6',	'TNFAIP8',	'TNFRSF9',	'TNFSF9',	'TNIP1',	'TNIP2',	'TRAF1',	
                       'TRIB1',	'TRIP10',	'TSC22D1',	'TUBB2A',	'VEGFA',	'YRDC',	'ZBTB10',	'ZC3H12A',	'ZFP36')
INFLAMMATORY_H <- c(	'ABCA1',	'ABI1',	'ACVR1B',	'ACVR2A',	'ADM',	'ADORA2B',	'ADRM1',	'AHR',	'APLNR',	'AQP9',	'ATP2A2',	'ATP2B1',	'ATP2C1',	'AXL',	
                     'BDKRB1',	'BEST1',	'BST2',	'BTG2',	'C3AR1',	'C5AR1',	'CALCRL',	'CCL17',	'CCL2',	'CCL20',	'CCL22',	'CCL24',	'CCL5',	'CCL7',	'CCR7',	
                     'CCRL2',	'CD14',	'CD40',	'CD48',	'CD55',	'CD69',	'CD70',	'CD82',	'CDKN1A',	'CHST2',	'CLEC5A',	'CMKLR1',	'CSF1',	'CSF3',	'CSF3R',	'CX3CL1',	
                     'CXCL10',	'CXCL11',	'CXCL6',	'CXCL9',	'CXCR6',	'CYBB',	'DCBLD2',	'EBI3',	'EDN1',	'EIF2AK2',	'EMP3',	'EMR1',	'EREG',	'F3',	'FFAR2',	'FPR1',	
                     'FZD5',	'GABBR1',	'GCH1',	'GNA15',	'GNAI3',	'GP1BA',	'GPC3',	'GPR132',	'GPR183',	'HAS2',	'HBEGF',	'HIF1A',	'HPN',	'HRH1',	'ICAM1',	'ICAM4',	
                     'ICOSLG',	'IFITM1',	'IFNAR1',	'IFNGR2',	'IL10',	'IL10RA',	'IL12B',	'IL15',	'IL15RA',	'IL18',	'IL18R1',	'IL18RAP',	'IL1A',	'IL1B',	'IL1R1',	'IL2RB',	
                     'IL4R',	'IL6',	'IL7R',	'IL8',	'INHBA',	'IRAK2',	'IRF1',	'IRF7',	'ITGA5',	'ITGB3',	'ITGB8',	'KCNA3',	'KCNJ2',	'KCNMB2',	'KIF1B',	'KLF6',	
                     'LAMP3',	'LCK',	'LCP2',	'LDLR',	'LIF',	'LPAR1',	'LTA',	'LY6E',	'LYN',	'MARCO',	'MEFV',	'MEP1A',	'MET',	'MMP14',	'MSR1',	'MXD1',	'MYC',	
                     'NAMPT',	'NDP',	'NFKB1',	'NFKBIA',	'NLRP3',	'NMI',	'NMUR1',	'NOD2',	'NPFFR2',	'OLR1',	'OPRK1',	'OSM',	'OSMR',	'P2RX4',	'P2RX7',	'P2RY2',	
                     'PCDH7',	'PDE4B',	'PDPN',	'PIK3R5',	'PLAUR',	'PROK2',	'PSEN1',	'PTAFR',	'PTGER2',	'PTGER4',	'PTGIR',	'PTPRE',	'PVR',	'RAF1',	'RASGRP1',	
                     'RELA',	'RGS1',	'RGS16',	'RHOG',	'RIPK2',	'RNF144B',	'ROS1',	'RTP4',	'SCARF1',	'SCN1B',	'SELE',	'SELL',	'SELS',	'SEMA4D',	'SERPINE1',	'SGMS2',	
                     'SLAMF1',	'SLC11A2',	'SLC1A2',	'SLC28A2',	'SLC31A1',	'SLC31A2',	'SLC4A4',	'SLC7A1',	'SLC7A2',	'SPHK1',	'SRI',	'STAB1',	'TACR1',	'TACR3',
                     'TAPBP',	'TIMP1',	'TLR1',	'TLR2',	'TLR3',	'TNFAIP6',	'TNFRSF1B',	'TNFRSF9',	'TNFSF10',	'TNFSF15',	'TNFSF9',	'TPBG',	'VIP')
HYPOXIA_H <- c(	'ADM',	'ADORA2B',	'AK4',	'AKAP12',	'ALDOA',	'ALDOB',	'ALDOC',	'AMPD3',	'ANGPTL4',	'ANKZF1',	'ANXA2',	'ATF3',	'ATP7A',	'B3GALT6',	'B4GALNT2',	
                'BCAN',	'BCL2',	'BGN',	'BHLHE40',	'BNIP3L',	'BRS3',	'BTG1',	'CA12',	'CASP6',	'CAV1',	'CCNG2',	'CCRN4L',	'CDKN1A',	'CDKN1B',	'CDKN1C',	'CHST2',	
                'CHST3',	'CITED2',	'COL5A1',	'CP',	'CSRP2',	'CTGF',	'CXCR4',	'CXCR7',	'CYR61',	'DCN',	'DDIT3',	'DDIT4',	'DPYSL4',	'DTNA',	'DUSP1',	'EDN2',	
                'EFNA1',	'EFNA3',	'EGFR',	'ENO1',	'ENO2',	'ENO3',	'ERO1L',	'ERRFI1',	'ETS1',	'EXT1',	'F3',	'FAM162A',	'FBP1',	'FOS',	'FOSL2',	'FOXO3',	'GAA',
                'GALK1',	'GAPDH',	'GAPDHS',	'GBE1',	'GCK',	'GCNT2',	'GLRX',	'GPC1',	'GPC3',	'GPC4',	'GPI',	'GRHPR',	'GYS1',	'HAS1',	'HDLBP',	'HEXA',	'HK1',	
                'HK2',	'HMOX1',	'HOXB9',	'HS3ST1',	'HSPA5',	'IDS',	'IER3',	'IGFBP1',	'IGFBP3',	'IL6',	'ILVBL',	'INHA',	'IRS2',	'ISG20',	'JMJD6',	'JUN',	
                'KDELR3',	'KDM3A',	'KIF5A',	'KLF6',	'KLF7',	'KLHL24',	'LALBA',	'LARGE',	'LDHA',	'LDHC',	'LOX',	'LXN',	'MAFF',	'MAP3K1',	'MIF',	'MT1E',	'MT2A',	
                'MXI1',	'MYH9',	'NAGK',	'NCAN',	'NDRG1',	'NDST1',	'NDST2',	'NEDD4L',	'NFIL3',	'NR3C1',	'P4HA1',	'P4HA2',	'PAM',	'PCK1',	'PDGFB',	'PDK1',	
                'PDK3',	'PFKFB3',	'PFKL',	'PFKP',	'PGAM2',	'PGF',	'PGK1',	'PGM1',	'PGM2',	'PHKG1',	'PIM1',	'PKLR',	'PKP1',	'PLAC8',	'PLAUR',	'PLIN2',	'PNRC1',	
                'PPARGC1A',	'PPFIA4',	'PPP1R15A',	'PPP1R3C',	'PRDX5',	'PRKCA',	'PRKCDBP',	'PTRF',	'PYGM',	'RBPJ',	'RORA',	'RRAGD',	'S100A4',	'SAP30',	'SCARB1',	
                'SDC2',	'SDC3',	'SDC4',	'SELENBP1',	'SERPINE1',	'SIAH2',	'SLC25A1',	'SLC2A1',	'SLC2A3',	'SLC2A5',	'SLC37A4',	'SLC6A6',	'SRPX',	'STBD1',	'STC1',	
                'STC2',	'SULT2B1',	'TES',	'TGFB3',	'TGFBI',	'TGM2',	'TIPARP',	'TKTL1',	'TMEM45A',	'TNFAIP3',	'TPBG',	'TPD52',	'TPI1',	'TPST2',	'UGP2',	
                'VEGFA',	'VHL',	'VLDLR',	'WISP2',	'WSB1',	'XPNPEP1',	'ZFP36',	'ZNF292')
EMT_H <- c(	'ABI3BP',	'ACTA2',	'ADAM12',	'ANPEP',	'APLP1',	'AREG',	'BASP1',	'BDNF',	'BGN',	'BMP1',	'CADM1',	'CALD1',	'CALU',	'CAP2',	'CAPG',	'CD44',	'CD59',	
            'CDH11',	'CDH2',	'CDH6',	'COL11A1',	'COL12A1',	'COL16A1',	'COL1A1',	'COL1A2',	'COL3A1',	'COL4A1',	'COL4A2',	'COL5A1',	'COL5A2',	'COL5A3',	'COL6A2',	'COL6A3',	
            'COL7A1',	'COL8A2',	'COMP',	'COPA',	'CRLF1',	'CTGF',	'CTHRC1',	'CXCL1',	'CXCL12',	'CXCL6',	'CYR61',	'DAB2',	'DCN',	'DKK1',	'DPYSL3',	'DST',	'ECM1',	'ECM2',	
            'EDIL3',	'EFEMP2',	'ELN',	'EMP3',	'ENO2',	'FAP',	'FAS',	'FBLN1',	'FBLN2',	'FBLN5',	'FBN1',	'FBN2',	'FERMT2',	'FGF2',	'FLNA',	'FMOD',	'FN1',	'FOXC2',	
            'FSTL1',	'FSTL3',	'FUCA1',	'FZD8',	'GADD45A',	'GADD45B',	'GAS1',	'GEM',	'GJA1',	'GLIPR1',	'GLT25D1',	'GPC1',	'GPX7',	'GREM1',	'HTRA1',	'ID2',	'IGFBP2',	
            'IGFBP3',	'IGFBP4',	'IL15',	'IL32',	'IL6',	'IL8',	'INHBA',	'ITGA2',	'ITGA5',	'ITGAV',	'ITGB1',	'ITGB3',	'ITGB5',	'JUN',	'LAMA1',	'LAMA2',	'LAMA3',	
            'LAMC1',	'LAMC2',	'LEPRE1',	'LGALS1',	'LOX',	'LOXL1',	'LOXL2',	'LRP1',	'LRRC15',	'LUM',	'MAGEE1',	'MATN2',	'MATN3',	'MCM7',	'MEST',	'MFAP5',	'MGP',	
            'MMP1',	'MMP14',	'MMP2',	'MMP3',	'MSX1',	'MXRA5',	'MYL9',	'MYLK',	'NID2',	'NNMT',	'NOTCH2',	'NT5E',	'NTM',	'OXTR',	'PCOLCE',	'PCOLCE2',	'PDGFRB',	'PDLIM4',	
            'PFN2',	'PLAUR',	'PLOD1',	'PLOD2',	'PLOD3',	'PMEPA1',	'PMP22',	'POSTN',	'PPIB',	'PRRX1',	'PRSS2',	'PTHLH',	'PTX3',	'PVR',	'QSOX1',	'RGS4',	'RHOB',	'SAT1',
            'SCG2',	'SDC1',	'SDC4',	'SERPINE1',	'SERPINE2',	'SERPINH1',	'SFRP1',	'SFRP4',	'SGCB',	'SGCD',	'SGCG',	'SLC6A8',	'SLIT2',	'SLIT3',	'SNAI2',	'SNTB1',	'SPARC',
            'SPOCK1',	'SPP1',	'TAGLN',	'TFPI2',	'TGFB1',	'TGFBI',	'TGFBR3',	'TGM2',	'THBS1',	'THBS2',	'THY1',	'TIMP1',	'TIMP3',	'TNC',	'TNFAIP3',	'TNFRSF11B',
            'TNFRSF12A',	'TPM1',	'TPM2',	'TPM4',	'VCAM1',	'VCAN',	'VEGFA',	'VEGFC',	'VIM',	'WIPF1',	'WNT5A')
SKIN_DEV <- c('CDH3',	'ABCB6',	'CDKN1A',	'ZMPSTE24',	'CDSN',	'YAP1',	'FST',	'TRIM16',	'TXNIP',	'PTGES3',	'TRAF3IP2',	'HPSE',	'EDAR',	'SLC27A4',	'SPINK5',	'SOX21',
                  'KRT71',	'MYSM1',	'KRT74',	'ASCL4',	'LCE7A',	'ACER1',	'KDF1',	'COL1A1',	'COL1A2',	'COL3A1',	'COL5A1',	'COL5A2',	'REG3G',	'COMP',	'CD109',	'CLDN4',	
                  'KRT72',	'KRT80',	'KRT25',	'DSG4',	'APCDD1',	'CSTA',	'CTNNB1',	'ASPRV1',	'CYP27B1',	'KRT28',	'SPRR4',	'DACT2',	'DHCR24',	'DLX3',	'DNASE1L2',	'JAG1',
                  'DSP',	'EDA',	'EGFR',	'KRT78',	'EPHA2',	'LCE4A',	'CERS3',	'ERCC2',	'EREG',	'ETV4',	'EVPL',	'EXT1',	'EXTL3',	'EZH2',	'FGF7',	'FGF10',	'FGFR2',
                  'DKK1',	'FOXC1',	'PALLD',	'FOXE1',	'EXPH5',	'FLG',	'FLNB',	'KAZN',	'PUM2',	'FOSL2',	'CASP14',	'OPN3',	'ALOX12',	'ALOX12B',	'ALOX15B',	'LCE5A',
                  'KLK5',	'POU2F3',	'SOSTDC1',	'CLIC4',	'SLITRK5',	'TRPC4AP',	'ABCA12',	'LCE2B',	'GATA6',	'VPS33B',	'GBA1',	'LATS2',	'GJB3',	'DKK4',	'INTU',	'AHDC1',
                  'GLI2',	'SFN',	'DLL1',	'KRT6C',	'GRHL1',	'SLC39A2',	'ANXA1',	'HDAC1',	'HDAC2',	'KRT73',	'HOXA7',	'HOXC13',	'KRTAP6-1',	'KRTAP6-2',	'KRTAP6-3',	
                  'KRT79',	'ZDHHC21',	'LIPM',	'KRT27',	'FOXI3',	'IGFBP5',	'RBPJ',	'LCE1A',	'LCE1B',	'LCE1C',	'LCE1D',	'LCE1E',	'LCE1F',	'LCE2A',	'LCE2C',	'LCE2D',
                  'LCE3A',	'LCE3B',	'LCE3C',	'LCE3E',	'IL1A',	'AQP3',	'IL17A',	'IL18',	'INHBA',	'ITGA6',	'IRF6',	'ITGA2',	'ITGA3',	'ITGB4',	'ITGB6',	'IVL',	'JUP',
                  'KRT77',	'CYSRT1',	'KRT1',	'KRT2',	'KRT3',	'KRT4',	'KRT5',	'KRT6A',	'KRT6B',	'KRT7',	'KRT9',	'KRT10',	'KRT14',	'KRT16',	'KRT17',	'HRNR',	'FLG2',	'KRT81',
                  'KRT82',	'KRT83',	'KRT84',	'KRT85',	'KRT86',	'TMPRSS11F',	'LAMA5',	'STMN1',	'LORICRIN',	'LRP4',	'LTB',	'MIR125B1',	'SMAD4',	'MET',	'ASAH1',	'MSX2',	
                  'LCE6A',	'MYD88',	'NAGLU',	'NF1',	'NGFR',	'NME2',	'NOTCH1',	'NUMA1',	'OVOL1',	'COL5A3',	'IL20',	'REG3A',	'PAX6',	'NSDHL',	'GAL',	'KRT76',	'WNT16',	
                  'PPHLN1',	'PDGFA',	'PIAS4',	'LSR',	'SUFU',	'ATP8A2',	'PKD1',	'PLEC',	'ATP7A',	'ERRFI1',	'SOX18',	'POU3F1',	'MED1',	'NSUN2',	'PPL',	'PPP3CA',	'LGR4',	
                  'TNFRSF19',	'MACROH2A2',	'FERMT1',	'PRKCH',	'ASH1L',	'MAP2K1',	'CYP26B1',	'PSEN1',	'ARRDC3',	'WDR48',	'GRHL3',	'OVOL2',	'PLAAT4',	'ALOXE3',	'BCL2',	
                  'RELA',	'SAV1',	'ALX4',	'ROCK1',	'BCR',	'RYR1',	'S100A7',	'SRSF6',	'NFKBIZ',	'LIPK',	'LIPN',	'NOM1',	'SHH',	'ELOVL1',	'BCL11B',	'SMO',	'SNAI1',	'SOS1',	
                  'SOX9',	'SPRR1A',	'SPRR1B',	'SPRR2B',	'SPRR2D',	'SPRR2E',	'SPRR2F',	'SPRR2G',	'SPRR3',	'SRF',	'ST14',	'ZFP36L1',	'STK4',	'TFAP2B',	'TGFB2',	'TGM1',	'TGM3',	
                  'TCHH',	'TNF',	'TSG101',	'UGCG',	'VDR',	'WNT5A',	'WNT10B',	'ZFP36',	'FA2H',	'SLC39A7',	'ZBED2',	'FZD3',	'GRHL2',	'FRAS1',	'IFT74',	'FUZ',	'WNT10A',	
                  'SLC2A10',	'SGPP1',	'SHARPIN',	'NCOA3',	'FZD6',	'CASP3',	'PIP5K1A',	'PLA2G10',	'TMEM79',	'CNFN',	'FOXN1',	'LCE3D',	'AP3B1',	'LGR5',	'TP63',	'PTCH2',	
                  'AKR1C3',	'KRT36',	'TRADD',	'ADAM9',	'SCEL',	'CFLAR',	'HDAC3',	'LDB1',	'CLDN1',	'LDB2',	'ACVR1B',	'LATS1',	'KRT75',	'GORAB',	'KLF4',	'LHX2',	'FOXQ1',	
                  'ROCK2',	'ADAMTS2',	'MACROH2A1',	'MAFB')
ANDROGEN_RESP_H <- c(	'ABCC4',	'ABHD2',	'ACSL3',	'ACTN1',	'ADAMTS1',	'ADRM1',	'AKAP12',	'AKT1',	'ALDH1A3',	'ANKH',	'APPBP2',	'ARID5B',	'AZGP1',	'B2M',	'B4GALT1',
                      'BMPR1B',	'CAMKK2',	'CCND1',	'CCND3',	'CDC14B',	'CDK6',	'CENPN',	'DBI',	'DHCR24',	'DNAJB9',	'ELK4',	'ELL2',	'ELOVL5',	'FADS1',	'FKBP5',	'GNAI3',	
                      'GPD1L',	'GSR',	'GUCY1A3',	'H1F0',	'HERC3',	'HMGCR',	'HMGCS1',	'HOMER2',	'HPGD',	'HSD17B14',	'IDI1',	'INPP4B',	'INSIG1',	'IQGAP2',	'ITGAV',	'KLK2',	
                      'KLK3',	'KRT19',	'KRT8',	'LIFR',	'LMAN1',	'MAF',	'MAK',	'MAP7',	'MERTK',	'MYL12A',	'NCOA4',	'NDRG1',	'NGLY1',	'NKX3-1',	'PA2G4',	'PDLIM5',	'PGM3',	
                      'PIAS1',	'PMEPA1',	'PPAP2A',	'PTK2B',	'PTPN21',	'RAB4A',	'RPS6KA3',	'RRP12',	'SAT1',	'SCD',	'SEC24D',	'SEPP1',	'SGK1',	'SLC26A2',	'SLC38A2',	'SMS',
                      'SORD',	'SPCS3',	'SPDEF',	'SRF',	'SRP19',	'STEAP4',	'STK39',	'TARP',	'TMEM50A',	'TMPRSS2',	'TNFAIP8',	'TPD52',	'TSC22D1',	'UAP1',	'UBE2I',	'UBE2J1',
                      'VAPA',	'XRCC5',	'XRCC6',	'ZBTB10',	'ZMIZ1',	'GRHL1',	'SLC39A2',	'ANXA1',	'HDAC1',	'HDAC2',	'KRT73',	'HOXA7',	'HOXC13',	'KRTAP6-1',	'KRTAP6-2',	
                      'KRTAP6-3',	'KRT79',	'ZDHHC21',	'LIPM',	'KRT27',	'FOXI3',	'IGFBP5',	'RBPJ',	'LCE1A',	'LCE1B',	'LCE1C',	'LCE1D',	'LCE1E',	'LCE1F',	'LCE2A',	'LCE2C',
                      'LCE2D',	'LCE3A',	'LCE3B',	'LCE3C',	'LCE3E',	'IL1A',	'AQP3',	'IL17A',	'IL18',	'INHBA',	'ITGA6',	'IRF6',	'ITGA2',	'ITGA3',	'ITGB4',	'ITGB6',	'IVL',	
                      'JUP',	'KRT77',	'CYSRT1',	'KRT1',	'KRT2',	'KRT3',	'KRT4',	'KRT5',	'KRT6A',	'KRT6B',	'KRT7',	'KRT9',	'KRT10',	'KRT14',	'KRT16',	'KRT17',	'HRNR',	'FLG2',	
                      'KRT81',	'KRT82',	'KRT83',	'KRT84',	'KRT85',	'KRT86',	'TMPRSS11F',	'LAMA5',	'STMN1',	'LORICRIN',	'LRP4',	'LTB',	'MIR125B1',	'SMAD4',	'MET',	'ASAH1',
                      'MSX2',	'LCE6A',	'MYD88',	'NAGLU',	'NF1',	'NGFR',	'NME2',	'NOTCH1',	'NUMA1',	'OVOL1',	'COL5A3',	'IL20',	'REG3A',	'PAX6',	'NSDHL',	'GAL',	'KRT76',	'WNT16',
                      'PPHLN1',	'PDGFA',	'PIAS4',	'LSR',	'SUFU',	'ATP8A2',	'PKD1',	'PLEC',	'ATP7A',	'ERRFI1',	'SOX18',	'POU3F1',	'MED1',	'NSUN2',	'PPL',	'PPP3CA',	'LGR4',	'TNFRSF19',
                      'MACROH2A2',	'FERMT1',	'PRKCH',	'ASH1L',	'MAP2K1',	'CYP26B1',	'PSEN1',	'ARRDC3',	'WDR48',	'GRHL3',	'OVOL2',	'PLAAT4',	'ALOXE3',	'BCL2',	'RELA',	'SAV1',	'ALX4',	
                      'ROCK1',	'BCR',	'RYR1',	'S100A7',	'SRSF6',	'NFKBIZ',	'LIPK',	'LIPN',	'NOM1',	'SHH',	'ELOVL1',	'BCL11B',	'SMO',	'SNAI1',	'SOS1',	'SOX9',	'SPRR1A',	'SPRR1B',	
                      'SPRR2B',	'SPRR2D',	'SPRR2E',	'SPRR2F',	'SPRR2G',	'SPRR3',	'SRF',	'ST14',	'ZFP36L1',	'STK4',	'TFAP2B',	'TGFB2',	'TGM1',	'TGM3',	'TCHH',	'TNF',	'TSG101',	'UGCG',
                      'VDR',	'WNT5A',	'WNT10B',	'ZFP36',	'FA2H',	'SLC39A7',	'ZBED2',	'FZD3',	'GRHL2',	'FRAS1',	'IFT74',	'FUZ',	'WNT10A',	'SLC2A10',	'SGPP1',	'SHARPIN',	'NCOA3',
                      'FZD6',	'CASP3',	'PIP5K1A',	'PLA2G10',	'TMEM79',	'CNFN',	'FOXN1',	'LCE3D',	'AP3B1',	'LGR5',	'TP63',	'PTCH2',	'AKR1C3',	'KRT36',	'TRADD',	'ADAM9',	'SCEL',	
                      'CFLAR',	'HDAC3',	'LDB1',	'CLDN1',	'LDB2',	'ACVR1B',	'LATS1',	'KRT75',	'GORAB',	'KLF4',	'LHX2',	'FOXQ1',	'ROCK2',	'ADAMTS2',	'MACROH2A1',	'MAFB')
# immune phenotypes
M1_score <-	c('NOS2',	'IL12A',	'FCGR1A',	'FCGR1B',	'CD80',	'CXCL9',	'CXCL10',	'CXCL11',
              'CD86',	'IL1A',	'IL1B',	'IL6',	'TNF',	'HLA-DRA',	'CCL5',	'IRF5',	'IRF1',	'CD40',	'IDO1',	'KYNU',	'CCR7')																							
M2_score <-	c('ARG1',	'IL10',	'FCGR2A',	'CD163',	'FCER2',	'CD200R1',	'PDCD1LG2',	'CD274',	'MARCO',	'CSF1R',	'MRC1',	
              'IL1RN',	'IL1R2',	'IL4R',	'CCL4',	'CCL13',	'CCL20',	'CCL17',	'CCL18',	'CCL22',	'CCL24',	'LYVE1',	
              'VEGFA',	'VEGFB',	'VEGFC',	'VEGFD',	'EGF',	'CTSA',	'CTSB',	'CSTC',	'CTSD',	'TGFB1',	'TGFB2',	'TGFB3',	'MMP14',	'MMP19',	
              'MMP9',	'CLEC7A',	'WNT7B',	'FASLG',	'TNFSF12',	'TNFSF8',	'CD276',	'VTCN1',	'MSR1',	'FN1',	'IRF4')
Cytotoxicity <-	c('GZMA',	'GZMB',	'GZMM',	'GZMK',	'GZMH',	'PRF1')
Exhaustion <- 	c('TOX',	'NR4A1',	'NR4A2',	'NR4A3',	'PDCD1',	'CTLA4',	'TIGIT',	'HAVCR2',	'LAG3',	'BTLA',	
                 'CD244',	'TBX21',	'EOMES',	'BLIMP1',	'NFAT',	'BATF',	'VHL',	'FOXO',	'FOXP1',	'ENTPD1',	
                 'CXCL13',	'MYO1E',	'PRDM1',	'WARS',	'CREM',	'PRDM1',	'LAYN',	'PHLDA1',	'SNAP47',	'CD38',	
                 'MYO7A',	'PARK7')	

#-------------------------------------------------------------------------------
# Add Module Score - Fig 7C-E

# Macrophages
Mac_seuratObj <- AddModuleScore(Mac_seuratObj,
                                features = list(M1_score),
                                name="M1_score")
VlnPlot(Mac_seuratObj, features = "M1_score1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) 

Mac_seuratObj <- AddModuleScore(Mac_seuratObj,
                                features = list(M2_score),
                                name="M2_score")
VlnPlot(Mac_seuratObj, features = "M2_score1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)  

# CTLs
CTL_seuratObj <- AddModuleScore(CTL_seuratObj,
                                features = list(Cytotoxicity),
                                name="Cytotoxicity")
VlnPlot(CTL_seuratObj, features = "Cytotoxicity1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) 

CTL_seuratObj <- AddModuleScore(CTL_seuratObj,
                                features = list(Exhaustion),
                                name="Exhaustion")
VlnPlot(CTL_seuratObj, features = "Exhaustion1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) 

# EXPRESSION OF MSIGDB GSEA PATHWAYS OF INTEREST IN EPITHLIAL CELLS 
Tumor_Epi_seuratObj <- AddModuleScore(Tumor_Epi_seuratObj,
                                features = list(TGFB_signaling_H),
                                name="TGFB_signaling_H")
VlnPlot(Tumor_Epi_seuratObj, features = "TGFB_signaling_H1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)  

Tumor_Epi_seuratObj <- AddModuleScore(Tumor_Epi_seuratObj,
                                      features = list(INFLAMMATORY_H),
                                      name="INFLAMMATORY_H")
VlnPlot(Tumor_Epi_seuratObj, features = "INFLAMMATORY_H1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) 

Tumor_Epi_seuratObj <- AddModuleScore(Tumor_Epi_seuratObj,
                                      features = list(HYPOXIA_H),
                                      name="HYPOXIA_H")
VlnPlot(Tumor_Epi_seuratObj, features = "HYPOXIA_H1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) 

Tumor_Epi_seuratObj <- AddModuleScore(Tumor_Epi_seuratObj,
                                      features = list(TNFA_signaling_H),
                                      name="TNFA_signaling_H")
VlnPlot(Tumor_Epi_seuratObj, features = "TNFA_signaling_H1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) 

Tumor_Epi_seuratObj <- AddModuleScore(Tumor_Epi_seuratObj,
                                      features = list(EMT_H),
                                      name="EMT_H")
VlnPlot(Tumor_Epi_seuratObj, features = "EMT_H1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) 

Tumor_Epi_seuratObj <- AddModuleScore(Tumor_Epi_seuratObj,
                                      features = list(SKIN_DEV),
                                      name="SKIN_DEV")
VlnPlot(Tumor_Epi_seuratObj, features = "SKIN_DEV1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) 

Tumor_Epi_seuratObj <- AddModuleScore(Tumor_Epi_seuratObj,
                                      features = list(ANDROGEN_RESP_H),
                                      name="ANDROGEN_RESP_H")
VlnPlot(Tumor_Epi_seuratObj, features = "ANDROGEN_RESP_H1", group.by = "FOXA1status", pt.size = 0) +stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)  
