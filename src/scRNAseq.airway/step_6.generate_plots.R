# step_6.generate_plots.R
# --- process scRNAseq data from hPSC-derived immune-airway organoids---
# step 6: Generate plots in the manuscript
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

# load cellchat object
cellchat.armm_arms <- readRDS(file.path(infodir, 'cellchat.armm_arms.rds'))

# Figure S3B:
# UMAP of hPSC-derived immune-airway organoids
my.cluster.color <- brewer.pal(12,'Paired')[1:7]
names(my.cluster.color) <- 0:6

g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix="Cluster", 
               tcolor=my.cluster.color, 
               tlabel=TRUE, tsplit=FALSE, tptsize=0.6) 
ggsave(file.path(figdir, "Fig.S2B.png"), plot=g, width=7.5, height=6.5, dpi=300)

# Figure S3C:
# UMAP of marker gene for each cluster of hPSC-derived immune-airway organoids.
g <- FeaturePlot(panc, features=c("KRT5","TUBA1A","CD14","KRT18","CoV2-N","COL1A1","XIST"), ncol=7)
ggsave(file.path(figdir, "Fig.S2C.png"), width=28, height=4.5, dpi=300)

# Figure S3D:
# Individual UMAP of immune-airway organoids exposed to mock (ARM+M) or SARS-CoV-2 (MOI=0.1, ARM+S), 
# and airway organoids co-cultured with 293T cells exposed to SARS-CoV-2 (MOI=0.1, ART+S).
g <- myDimPlot3(tobj=panc, treduct="umap", tgroup_by="ident", tgroup_order=0:6, tsuffix="Cluster", 
                tcolor=my.cluster.color, tsplit_by='Name', 
                tsplit_order=c('ARMM','ARMS','ARTS'), 
                tlabel=T, tncol=3, tptsize=0.4, tlbsize=5) + theme(legend.position="none")
ggsave(file.path(figdir, "Fig.S2D.png"), plot=g, width=12.5, height=4.5, dpi=300)

# Figure S3E:
# Dot plot analysis of proinflammatory macrophage associated genes in macrophage cluster of ARM+M or ARM+S conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c("IL1B","IL6","CD80"), 
                 cluster.id=2, name.order=c('ARMM','ARMS'))
ggsave(file.path(figdir, "Fig.S3E.png"), plot=g, width=4.5, height=3, dpi=300)

# Figure 3H:
# Dot plot analysis of SASP associate genes in Ciliated cell cluster of ARM+S and ARM+T conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c(c("ANAPC16","ANAPC7","ATM","BMI1","CCNE2","CDK4",
                                            "CDKN2B","HMGA2","MAPK14","MAPKAPK2","MAPKAPK3","NFKB1",
                                            "RBBP7","TERF1","TFDP1","TFDP2","TNRC6A","TP53")), 
                 cluster.id=1, name.order=c('ARTS','ARMS'))
ggsave(file.path(figdir, "Fig.S3H.png"), plot=g, width=4.5, height=6.5, dpi=300)

# Figure 3I:
# Dot plot analysis of senescence associate genes in Ciliated cell cluster of macrophage-airway organoids 
# of ARM+S and ARM+T conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c("TNFRSF1B","IL1A","CCL20","IL1B","IL6","XCL1",
                                          "VCAM1","CXCL5","PF4","CXCL1","SELL"), 
                 cluster.id=1, name.order=c('ARTS','ARMS'))
ggsave(file.path(figdir, "Fig.S3H.png"), plot=g, width=4.5, height=5, dpi=300)

# Figure S6B:
# Dot plot showing the differential signaling from macrophages to ciliated cells in 
# hPSC-derived ARM+M and ARM+S conditions.
g <- netVisual_bubble(cellchat.armm_arms, sources.use=c(5), targets.use = c(3), 
                      comparison = c(1,2), max.dataset = 2, 
                      title.name = "Increased signaling", 
                      angle.x = 45, remove.isolate = T, line.on=F)
ggsave(file.path(figdir, "Fig.S6B.png"), width=3.5, height=9, dpi=300)
