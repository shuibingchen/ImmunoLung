# step_4.cell_clustering.R
# --- process scRNAseq data from hPSC-derived immune-alveolar organoids ---
# step 4: cluster cells and infer cell types
# Author: Tuo Zhang
# Date: 12/1/2024
# 

library(Seurat)
library(scater)
library(scran)
library(batchelor)
library(tidyverse)
library(scuttle)
library(future)
library(RColorBrewer)
library(pheatmap)
library(R.utils)
library(magrittr)
library(patchwork)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# Set a random seed
set.seed(98)

# set parallelization in Seurat
plan("multicore", workers=6)
options(future.globals.maxSize=10*1024^3)

# Cluster the cells
my.dims <- 1:50
panc %<>% FindNeighbors(reduction="mnn", dims=my.dims)
panc %<>% FindClusters(resolution=0.8, verbose=T)

# Run UMAP dimensionality reduction
panc %<>% RunUMAP(dims=my.dims, reduction="mnn", n.components=2, seed.use=42, 
                  n.neighbors=30, n.epochs=500)

# set cell identity
panc %<>% SetIdent(value="RNA_snn_res.0.8")

# merge clusters
merge.clust <- c(2,0,0,1,0,1,0,3,1,0,4,1,2,5,3,6,3,5,3,1)
names(merge.clust) <- 0:19

final.clust <- as.vector(merge.clust[as.vector(Idents(panc))])
names(final.clust) <- names(Idents(panc))
final.clust <- factor(final.clust, levels=0:6)

# add final clusters to meta data
panc[["merged.cluster"]] <- final.clust

# set final clusters
Idents(panc) <- "merged.cluster"

# identify marker genes for each cluster
for (k in 0:6){
  tmarkers <- FindMarkers(object=panc, ident.1=k, min.pct=0.25, test.use="wilcox", assay="RNA")
  # write to file
  write.table(as.data.frame(tmarkers), file=file.path(infodir, paste0("C", k, ".markers.txt")), quote=FALSE, na="", sep="\t", col.names=NA)
}

# save Seurat object
saveRDS(panc, file.path(infodir, 'panc.rds'))
