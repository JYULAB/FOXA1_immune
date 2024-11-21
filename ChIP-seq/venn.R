BiocManager::install("VennDiagram",force = T)
library(VennDiagram)

P_FOXA1_peaks <- read.csv(file = "Luly_FOXA1/SKO_only_FOXA1_venn", header = T, sep = '')
PF_FOXA1_peaks <- read.csv(file = "Luly_FOXA1/DKO_only_FOXA1_venn", header = T, sep = '')

FOXA1_over <- venn.diagram(
  x = list(P_FOXA1_peaks$X.name,PF_FOXA1_peaks$X.name),
  category.names = c("SKO","DKO"),
  margin = 0.15,
  filename = NULL,
  output=T,
  imagetype="png" ,
  height = 1100 , 
  width = 1100 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#4994c4","#ea5414"),
  cex = 0.5,
  fontfamily = "sans",
  ext.dist = 0.1,
  ext.length = 0.8,
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(-90, 90),
  cat.dist = c(0.25, 0.25),
  cat.fontfamily = "sans",
  cat.col = c("#000000")
)
pdf(file="FOXA1_over_spike_in.pdf")
grid.draw(FOXA1_over)
dev.off()