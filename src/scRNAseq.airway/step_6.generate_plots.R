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

# Figure S4B:
# UMAP of hPSC-derived immune-airway organoids
my.cluster.color <- brewer.pal(12,'Paired')[1:7]
names(my.cluster.color) <- 0:6

g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix="Cluster", 
               tcolor=my.cluster.color, 
               tlabel=TRUE, tsplit=FALSE, tptsize=0.6) 
ggsave(file.path(figdir, "Fig.S4B.png"), plot=g, width=7.5, height=6.5, dpi=300)

# Figure S4C:
# UMAP of marker gene for each cluster of hPSC-derived immune-airway organoids.
g <- FeaturePlot(panc, features=c("KRT5","TUBA1A","CD14","KRT18","CoV2-N","COL1A1","XIST"), ncol=7)
ggsave(file.path(figdir, "Fig.S4C.png"), width=28, height=4.5, dpi=300)

# Figure S4D:
# Correlation analysis of genes with cell fates in hPSC-derived immuno-airway organoids and adult human lung cells. 
#   Load marker genes for our data
our.markers <- list()

for (k in 0:6){
  tmarker <- read.table(file.path(infodir, paste0("C", k, ".markers.txt")), 
                        header=T, check.names=F, stringsAsFactors=F, row.names=1, sep='\t')
  # select positive markers
  tmarker <- subset(tmarker, p_val_adj < 0.05 & avg_log2FC > 0)
  # add to list
  our.markers[[paste0('C',k)]] <- rownames(tmarker)
}

#   Load marker genes for reference data
ref.df <- read.table(file.path(sourcedir, "ref.markers.txt.gz"), 
                     header=T, check.names=F, stringsAsFactors=F, sep='\t')

ref.markers <- list()

for (c in unique(ref.df$celltype)){
  ref.markers[[c]] <- as.vector(ref.df %>% dplyr::filter(celltype==c) %>% pull(gene))
}

plotFracOverlap(ref.markers[c('Ciliated','Proximal Basal','Proliferating basal','Myofibroblast')], 
                out.markers[paste0('C', c(0,1,3,5))], 
                file.path(figdir, "Fig.S4D.png"), width=4, height=4)

# Figure S4E:
# UMAP of additional marker genes for each cluster of immuno-airway organoids.
g <- FeaturePlot(panc, features=c("FOXJ1","DNAH5","PIFO"))
ggsave(file.path(figdir, "Fig.S4E.png"), plot=g, width=12.5, height=4.5, dpi=300)

# Figure S5A:
# Individual UMAP of immune-airway organoids exposed to mock (ARM+M) or SARS-CoV-2 (MOI=0.1, ARM+S), 
# and airway organoids co-cultured with 293T cells exposed to SARS-CoV-2 (MOI=0.1, ART+S).
g <- myDimPlot3(tobj=panc, treduct="umap", tgroup_by="ident", tgroup_order=0:6, tsuffix="Cluster", 
                tcolor=my.cluster.color, tsplit_by='Name', 
                tsplit_order=c('ARMM','ARMS','ARTS'), 
                tlabel=T, tncol=3, tptsize=0.4, tlbsize=5) + theme(legend.position="none")
ggsave(file.path(figdir, "Fig.S5A.png"), plot=g, width=12.5, height=4.5, dpi=300)

# Figure S5B:
# Dot plot analysis of proinflammatory macrophage associated genes in macrophage cluster of ARM+M or ARM+S conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c("IL1B","IL6","CD80"), 
                 cluster.id=2, name.order=c('ARMM','ARMS'))
ggsave(file.path(figdir, "Fig.S5B.png"), plot=g, width=4.5, height=3, dpi=300)

# Figure S5E:
# Dot plot analysis of SASP associate genes in Ciliated cell cluster of ARM+S and ARM+T conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c(c("ANAPC16","ANAPC7","ATM","BMI1","CCNE2","CDK4",
                                            "CDKN2B","HMGA2","MAPK14","MAPKAPK2","MAPKAPK3","NFKB1",
                                            "RBBP7","TERF1","TFDP1","TFDP2","TNRC6A","TP53")), 
                 cluster.id=1, name.order=c('ARTS','ARMS'))
ggsave(file.path(figdir, "Fig.S5E.png"), plot=g, width=4.5, height=6.5, dpi=300)

# Figure S5F:
# Dot plot analysis of senescence associate genes in Ciliated cell cluster of macrophage-airway organoids 
# of ARM+S and ARM+T conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c("TNFRSF1B","IL1A","CCL20","IL1B","IL6","XCL1",
                                          "VCAM1","CXCL5","PF4","CXCL1","SELL"), 
                 cluster.id=1, name.order=c('ARTS','ARMS'))
ggsave(file.path(figdir, "Fig.S5F.png"), plot=g, width=4.5, height=5, dpi=300)

# Figure S8B:
# Dot plot showing the differential signaling from macrophages to ciliated cells in 
# hPSC-derived ARM+M and ARM+S conditions.
g <- netVisual_bubble(cellchat.armm_arms, sources.use=c(5), targets.use = c(3), 
                      comparison = c(1,2), max.dataset = 2, 
                      title.name = "Increased signaling", 
                      angle.x = 45, remove.isolate = T, line.on=F)
ggsave(file.path(figdir, "Fig.S8B.png"), width=3.5, height=9, dpi=300)

