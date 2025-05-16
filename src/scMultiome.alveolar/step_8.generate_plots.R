# step_8.generate_plots.R
# --- process multiome (RNA+ATAC) data from hPSC-derived immune-alveolar organoids---
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
cellchat.avum_avpm <- readRDS(file.path(infodir, 'cellchat.avum_avpm.rds'))

# set color
my.cluster.color <- brewer.pal(12,'Paired')[1:8]
names(my.cluster.color) <- 0:7

# Figure 3A: 
# UMAP of hPSC-derived immune-alveolar organoids containing unstimulated (AVUM) and 
# proinflammatory macrophages (AVPM).
g <- myDimPlot(tobj=panc, treduct="wnn.umap", tcate="ident", torder=0:7, 
                    tcolor=my.cluster.color, tsuffix="Cluster", tlabel=TRUE, 
                    tsplit=FALSE, tptsize=0.6)
ggsave(file.path(figdir, "Fig.3A.png"), plot=g, width=8, height=6, dpi=300)

# Figure 3B:
# Dot plot displaying cell marker genes.
g <- DotPlot(panc, features=c("SFTPB","CLIC5","KRT8","CD163","COL1A1")) + coord_flip()
ggsave(file.path(figdir, "Fig.3B.png"), width=7.5, height=4, dpi=300)

# Figure 3C:
# Individual UMAP of snRNA-seq and snATAC-seq analysis of hPSC-derived AVUM and AVPM organoids.
p1 <- myDimPlot(tobj=panc, treduct="rna_umap", tcate="ident", torder=0:7, 
                tcolor=my.cluster.color, tsuffix="RNA", tlabel=TRUE, 
                tsplit=FALSE, tptsize=0.6) + theme(legend.position="bottom")

p2 <- myDimPlot(tobj=panc, treduct="atac_integrated_umap", tcate="ident", torder=0:7,
                tcolor=my.cluster.color, tsuffix="ATAC", tlabel=TRUE, 
                tsplit=FALSE, tptsize=0.6) + theme(legend.position="bottom")

g <- p1 | p2
ggsave(file.path(figdir, "Fig.3C.png"), plot=g, width=12.5, height=6.5, dpi=300)

# Figure 3D: 
# Dot plot analysis of proinflammatory macrophage-associated genes in macrophage cluster of 
# hPSC-derived AVUM and AVPM organoids. 
g <- my.dot.plot(seurat.obj=panc, genes=c("IL1B","IL6","CD80","IDO1"), 
                 cluster.id=3, name.order=c('AVUM','AVPM'))
ggsave(file.path(figdir, "Fig.3D.png"), plot=g, width=4.5, height=3.5)

# Figure 3E:
# Dot plot analysis of SASP associate genes in AT2 cell cluster of hPSC-derived AVUM and AVPM organoids. 
g <- my.dot.plot(seurat.obj=panc, genes=c("CSF2RB","CCL17","CSF3","IL1B","IL6","CXCL1","CCL27",
                                          "CXCL5","IGFBP3","CXCL2"), 
                 cluster.id=0, name.order=c('AVUM','AVPM'))
ggsave(file.path(figdir, "Fig.3E.png"), plot=g, width=5, height=4)

# Figure 3F:
# Dot plot analysis of senescence associate genes in AT2 cell cluster of hPSC-derived AVUM and AVPM organoids.
g <- my.dot.plot(seurat.obj=panc, genes=c("ANAPC1","ANAPC16","CDK4","CDK6","ETS1","MAPK10","MAPKAPK2",
                                          "MAPKAPK3","NFKB1","UBB","TP53","CXCL8","SIRT1"), 
                 cluster.id=0, name.order=c('AVUM','AVPM'))
ggsave(file.path(figdir, "Fig.3F.png"), plot=g, width=5.5, height=4.5)

# Figure 4A:
# Dot plot showed the upregulated signals from macrophages to AT2 cells in 
# alveolar organoids containing unstimulated (AVUM) or proinflammatory macrophages (AVPM).
g <- netVisual_bubble(cellchat.avum_avpm, sources.use=c(7), targets.use = c(2), 
                      comparison = c(1,2), max.dataset = 2, 
                      title.name = "Increased signaling", 
                      angle.x = 45, remove.isolate = T, line.on=F)
ggsave(file.path(figdir, "Fig.4A.png"), width=3.5, height=4.5, dpi=300)

# Figure S6B
# Individual UMAP of snRNA-seq and snATAC-seq of AVUM and AVPM macrophages.
g <- myDimPlot3(tobj=panc, treduct="wnn.umap", tgroup_by="ident", tgroup_order=0:7, tsuffix="RNA_Cluster", 
                tcolor=my.cluster.color, tsplit_by='Name', 
                tsplit_order=c("AVUM","AVPM"), 
                tlabel=T, tncol=2, tptsize=0.4, tlbsize=4) + theme(legend.position="none")
ggsave(file.path(figdir, "Fig.S6B.png"), plot=g, width=8, height=4, dpi=300)

# Figure S6C
# Dot plot analysis of fibrosis associated genes in fibroblast cluster of AVUM and AVPM conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c("HIF1A","SKP2","MTOR","CXCR4","S1PR2","IL10",
                                          "CFLAR","ITGB1","ITGA2","ACTA2"), 
                 cluster.id=6, name.order=c('AVUM','AVPM'))
ggsave(file.path(figdir, "Fig.S6C.png"), plot=g, width=4.5, height=5)

