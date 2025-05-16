# step_7.generate_plots.R
# --- process scRNAseq data from hPSC-derived immune-alveolar organoids---
# step 7: Generate plots in the manuscript
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
cellchat.avmm_avms <- readRDS(file.path(infodir, 'cellchat.avmm_avms.rds'))

# load gene set enrichment results of
# Immune-alveolar organoids exposed to SARS-CoV-2 vs. 
# alveolar organoids co-cultured with 293T cells exposed to SARS-CoV-2
gse <- readRDS(file.path(infodir, "gse.DE.C0-AVMS.vs.C0-AVTS.mincells_10.wilcox.rds"))

# Figure 2A:
# UMAP of hPSC-derived immune-alveolar organoids analyzed by scRNA-seq.
my.cluster.color <- brewer.pal(12,'Paired')[1:8]
names(my.cluster.color) <- 0:7

g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix="Cluster", 
               tcolor=my.cluster.color, 
               tlabel=TRUE, tsplit=FALSE, tptsize=0.6) 
ggsave(file.path(figdir, "Fig.2A.png"), plot=g, width=7.5, height=6.5, dpi=300)

# Figure 2B:
# Dot plot displaying cell marker gene of each cluster of hPSC-derived immune-alveolar organoids.
g <- DotPlot(panc, features=c("SFTPB","KRT8","VEGFA","CD163","COL1A1","XIST")) + coord_flip()
ggsave(file.path(figdir, "Fig.2B.png"), width=6.5, height=4, dpi=300)

# Figure 2C:
# Individual UMAP of immune-alveolar organoids exposed to mock (AVM+M) or SARS-CoV-2 (MOI=0.1, AVM+S), 
# and alveolar organoids co-cultured with 293T cells exposed to SARS-CoV-2 (MOI=0.1, AVT+S).
g <- myDimPlot3(tobj=panc, treduct="umap", tgroup_by="ident", tgroup_order=0:7, tsuffix="Cluster", 
                tcolor=my.cluster.color, tsplit_by='Name', 
                tsplit_order=c('AVMM','AVMS','AVTS'), 
                tlabel=T, tncol=3, tptsize=0.4, tlbsize=5) + theme(legend.position="none")
ggsave(file.path(figdir, "Fig.2C.png"), plot=g, width=12.5, height=4.5, dpi=300)

# Figure 2D:
# Dot plot analysis of proinflammatory macrophage-associated genes in macrophage cluster of 
# AVM+M and AVM+S conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c("IL1B","IL6","IDO1","CD80"), 
            cluster.id=3, name.order=c('AVMM','AVMS'))
ggsave(file.path(figdir, "Fig.2D.png"), plot=g, width=4.5, height=4, dpi=300)

# Figure 2G: 
# Enrichment of cell death pathways in AT2 cell cluster of 
# immuno-alveolar organoids (AVM+S) or 293T co-cultured with alveolar organoids (AVT+S) exposed SARS-CoV-2 (MOI=0.1). 
# bar plot
gse.df <- as.data.frame(gse) 
description.order <- gse.df %>% arrange(NES) %>% pull(Description)
gse.df <- gse.df %>% 
  mutate(Sign=ifelse(enrichmentScore > 0, "activated", "suppressed")) %>% 
  mutate(Description=factor(Description, levels=description.order))
g <- ggplot(tdata, aes(x=Description, y=NES, fill=Sign))
g <- g + geom_bar(stat="identity")
g <- g + scale_fill_manual(values=c('activated'='#4575b4','suppressed'='#fc8d59'))
g <- g + geom_hline(yintercept = 0, linewidth=0.5)
g <- g + coord_flip()
g <- g + ylab("Normalized Enrichment Score")
g <- g + theme_bw()
g <- g + theme(legend.position="none", panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), axis.title.y=element_blank(),
               axis.text = element_text(color='black'),
               axis.title.x = element_text(color='black'))
ggsave(file.path(figdir, "Fig.2G.png"), width=5, height=4.5, dpi=300)

# Figure 2H:
# Gene Set Enrichment Analysis (GSEA) of senescence pathway in AT2 cell cluster of 
# AVM+S versus AVT+S condition.
pos <- which(gse$Description=="Senescence")
png(file.path(figdir, "Fig.2H.png"), units='in', width=8, height=5.5, res=300)
print(gseaplot2(gse, geneSetID = pos, title=gse$Description[pos]))
dev.off()

# Figure 2I:
# Dot plot analysis of senescence-associated secretory phenotype (SASP) associated genes in 
# AT2 cell cluster of AVM+S and AVT+S conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c("CSF2RB","CCL17","CSF3","IL1B","IL6","CXCL1","CCL27","CCL2","LEPR"), 
                 cluster.id=0, name.order=c('AVTS','AVMS'))
ggsave(file.path(figdir, "Fig.2I.png"), plot=g, width=5, height=4.5, dpi=300)

# Figure 2J:
# Dot plot analysis of senescence associate genes in AT2 cell cluster of 
# AVM+S and AVT+S conditions.
g <- my.dot.plot(seurat.obj=panc, 
                 genes=c("ANAPC1","ANAPC16","CDK4","CDK6","ETS1",
                         "MAPK10","MAPKAPK2","MAPKAPK3","NFKB1","UBB",
                         "CXCL8","TP53","SIRT1"), 
                 cluster.id=0, name.order=c('AVTS','AVMS'))
ggsave(file.path(figdir, "Fig.2J.png"), plot=g, width=5, height=5.5, dpi=300)

# Figure 4B:
# Dot plot showed the upregulated signals from macrophages to AT2 cells in 
# immune-alveolar organoids exposed to mock (AVM+M) versus SARS-CoV-2 (MOI=0.1, AVM+S).
g <- netVisual_bubble(cellchat.avmm_avms, sources.use=c(7), targets.use = c(2), 
                      comparison = c(1,2), max.dataset = 2, 
                      title.name = "Increased signaling", 
                      angle.x = 45, remove.isolate = T, line.on=F)
ggsave(file.path(figdir, "Fig.4B.png"), width=3.5, height=6.75, dpi=300)

# Figure S2D:
# UMAP of marker gene for each cluster of immune-alveolar organoids 
# exposed to mock (AVM+M) or SARS-CoV-2 (MOI=0.1, AVM+S), and alveolar organoids 
# co-cultured with 293T cells exposed to SARS-CoV-2 (MOI=0.1, AVT+S). 
g <- FeaturePlot(panc, features=c("SFTPB","KRT8","VEGFA","CD14","COL1A1","XIST"))
ggsave(file.path(figdir, "Fig.S2D.png"), width=12.5, height=8.5, dpi=300)

# Figure S2E:
# Correlation analysis of genes with cell fates in hPSC-derived immuno-alveolar organoids and adult human lung cells. 
#   Load marker genes for our data
our.markers <- list()

for (k in 0:7){
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

plotFracOverlap(ref.markers[c('Basal','Proximal Basal','Differentiating basal','Proliferating basal',
                              'Alveolar epithelial Type 1','Signaling alveolar epithelial Type 2',
                              'Myofibroblast')], 
                out.markers[paste0('C', c(0,1,2,4,5,6))], 
                file.path(figdir, "Fig.S2E.png"), width=4, height=5)

# Figure S2F:
# UMAP of additional marker genes for each cluster of immuno-alveolar organoids
g <- FeaturePlot(panc, features=c("NAPSA","KRT5","TP63","KRT14","ITGA6"))
ggsave(file.path(figdir, "Fig.S2F.png"), width=15.5, height=4.5, dpi=300)

# Figure S3A:
# Dot plot showing the detection of SARS-CoV2-M transcript hPSC-derived immuno-alveolar organoids 
# exposed to mock (AVM+M) or SARS-CoV-2 (MOI=0.25, AVM+S), and 
# alveolar organoids co-cultured with 293T cells exposed to SARS-CoV-2 (MOI=0.25, AVT+S).  
panc$Name <- factor(panc$Name, levels=c('AVMM','AVMS','AVTS'))
g <- VlnPlot(panc, features="CoV2-M", group.by="Name")
ggsave(file.path(figdir, "Fig.S3A.png"), width=3.5, height=2.25, dpi=300)

# Figure S3B:
# Dot plot analysis of SASP associated genes in AT2 cell cluster of AVM+S and AVM+M conditions. 
g <- my.dot.plot(seurat.obj=panc, genes=c("CXCL1","IL6","CSF3","CCL17","CSF2RB"),
            cluster.id=0, name.order=c('AVMM','AVMS'))
ggsave(file.path(figdir, "Fig.S3B.png"), plot=g, width=5, height=3.5, dpi=300)

# Figure S3C:
# Dot plot analysis of senescence associate genes in AT2 cell cluster of AVM+S and AVM+M conditions.
g <- my.dot.plot(seurat.obj=panc, genes=c("NFKB1","MAPKAPK2","MAPK10","ETS1","CDK4","CDK6","ANAPC16","ANAPC1"),
                 cluster.id=0, name.order=c('AVMM','AVMS'))
ggsave(file.path(figdir, "Fig.S3C.png"), plot=g, width=5, height=5, dpi=300)

# Figure S3J:
# Dot plot analysis of senescence associated genes in basal cells of AVM+S and AVM+M conditions.
g <- my.dot.plot(seurat.obj=panc, 
                 genes=c("SIRT1","MMP1","CXCL8","UBE2C","UBB","TP53","SUZ12","SP1","MAPKAPK2","HMGA2","EZH2",
                         "EED","CDK4","CDK2","CDC27","CCNE2","CBX2","ATM","ANAPC1"),
                 cluster.id=c(1,4,5), name.order=c('AVTS','AVMS'))
ggsave(file.path(figdir, "Fig.S3J.png"), plot=g, width=5, height=5.5, dpi=300)

# Figure S3K:
# Dot plot analysis of fibrosis associated genes in fibroblast cluster of AVM+S and AVT+S conditions.
g <- my.dot.plot(seurat.obj=panc, 
                 genes=c("HIF1A","SKP2","MTOR","CXCR4","S1PR2","IL10","CFLAR","ITGB1","ITGA2","ACTA2"), 
                 cluster.id=6, name.order=c('AVTS','AVMS'))
ggsave(file.path(figdir, "Fig.S2K.png"), width=4.5, height=5, dpi=300)

