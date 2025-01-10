# step_8.generate_plots.R
# --- process multiome (RNA+ATAC) data from hPSC-derived immune-airway organoids---
# step 8: Generate plots in the manuscript
# Author: Tuo Zhang
# Date: 12/1/2024
# 

library(Seurat)
library(CellChat)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(R.utils)
library(magrittr)
library(patchwork)
library(clusterProfiler)
library(enrichplot)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# load Seurat object
panc <- readRDS(file.path(infodir, 'panc.rds'))

DefaultAssay(panc) <- "RNA"

# load cellchat object
cellchat.arum_arpm <- readRDS(file.path(infodir, 'cellchat.arum_arpm.rds'))

# set color
my.cluster.color <- brewer.pal(12,'Paired')[1:5]
names(my.cluster.color) <- 0:4

# Figure S5B: 
# UMAP of ARUM and ARPM organoids.
g <- myDimPlot(tobj=panc, treduct="wnn.umap", tcate="ident", torder=0:4, 
               tcolor=my.cluster.color, tsuffix="Cluster", tlabel=TRUE, 
               tsplit=FALSE, tptsize=0.6)
ggsave(file.path(figdir, "Fig.S5B.png"), plot=g, width=7.5, height=6.5, dpi=300)

# Figure S5C:
# UMAP of marker gene for each cluster of hPSC-derived immune-airway organoids.
g <- FeaturePlot(panc, features=c("TUBA1A","KRT17","COL1A1","CD163"), 
                 reduction="wnn.umap", ncol=2)
ggsave(file.path(figdir, "Fig.S5C.png"), plot=g, width=8.5, height=8, dpi=300)

# Figure S5D:
# UMAP projections of single-cell RNA (left) and ATAC (right) data of hPSC-derived immune-airway organoids,
# with cells colored by cell types based on transcriptomic profiles. 
p1 <- myDimPlot(tobj=panc, treduct="rna_umap", tcate="ident", torder=0:4, 
                tcolor=my.cluster.color, tsuffix="RNA", tlabel=TRUE, 
                tsplit=FALSE, tptsize=0.6) + theme(legend.position="bottom")

p2 <- myDimPlot(tobj=panc, treduct="atac_integrated_umap", tcate="ident", torder=0:4,
                tcolor=my.cluster.color, tsuffix="ATAC", tlabel=TRUE, 
                tsplit=FALSE, tptsize=0.6) + theme(legend.position="bottom")

g <- p1 | p2
ggsave(file.path(figdir, "Fig.S5D.png"), plot=g, width=12.5, height=6.5, dpi=300)

# Figure S5E:
# Joint UMAP projections of RNA and ATAC data showing varied cell types under ARUM (left) and ARPM (right) conditions.
g <- myDimPlot3(tobj=panc, treduct="wnn.umap", tgroup_by="ident", tgroup_order=0:4, 
                tsuffix="RNA_Cluster", tcolor=my.cluster.color, 
                tsplit_by='Name', tsplit_order=c("ARUM","ARPM"), 
                tlabel=T, tncol=2, tptsize=0.4, tlbsize=4) + theme(legend.position="none")
ggsave(file.path(figdir, "Fig.S5E.png"), plot=g, width=8, height=4, dpi=300)

# Figure S5F:
# Dot plot showing the expressions of pro-inflammatory macrophage associated genes in 
# macrophage cell cluster of ARUM and ARPM organoids. 
g <- my.dot.plot(seurat.obj=panc, genes=c("IL1B","IL6","CD80"), 
                 cluster.id=3, name.order=c('ARUM','ARPM'))
ggsave(file.path(figdir, "Fig.S6F.png"), plot=g, width=4.5, height=3.5)

# Figure S5G:
# Dot plot showing the expressions of senescence associated genes in 
# ciliated cell cluster of ARUM and ARPM organoids.
g <- my.dot.plot(seurat.obj=panc, genes=c("ANAPC16","ANAPC7","ATM","BMI1","CCNE2","CDK4","CDKN2B",
                                          "HMGA2","MAPK14","MAPKAPK2","MAPKAPK3","NFKB1","RBBP7",
                                          "TERF1","TFDP1","TFDP2","TNRC6A","TP53"), 
                 cluster.id=0, name.order=c('ARUM','ARPM'))
ggsave(file.path(figdir, "Fig.S6G.png"), plot=g, width=4.5, height=6.5)

# Figure S5H:
# Dot plot showing the expressions of senescence associated genes in 
# ciliated cell cluster of ARUM and ARPM organoids. 
g <- my.dot.plot(seurat.obj=panc, genes=c("TNFRSF1B","IL1A","CCL20","IL1B","IL6","XCL1","VCAM1",
                                          "CXCL5","PF4","CXCL1","SELL"), 
                 cluster.id=0, name.order=c('ARUM','ARPM'))
ggsave(file.path(figdir, "Fig.S6H.png"), plot=g, width=4.5, height=5)

# Figure S6A:
# Dot plot showing the differential signaling from macrophages to ciliated cells in ARUM and ARPM organoids.
g <- netVisual_bubble(cellchat.arum_arpm, sources.use=c(5), targets.use = c(3), 
                      comparison = c(1,2), max.dataset = 2, 
                      title.name = "Increased signaling", 
                      angle.x = 45, remove.isolate = T, line.on=F)
ggsave(file.path(figdir, "Fig.S6A.png"), width=3.5, height=4, dpi=300)
