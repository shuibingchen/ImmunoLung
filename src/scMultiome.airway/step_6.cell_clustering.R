# step_6.cell_clustering.R
# --- process multiome (RNA+ATAC) data from hPSC-derived immune-airway organoids ---
# step 6: Perform cell clustering based on RNA data, merge RNA and ATAC data
# Author: Tuo Zhang
# Date: 12/1/2024
# 

library(Seurat)
library(Signac)
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
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# random seed
set.seed(98)

# load functions
setwd(workdir)
source("my_functions.R")

# set parallelization in Seurat
plan("multicore", workers=6)
options(future.globals.maxSize=10*1024^3)

# load RNA Seurat object
panc.rna <- readRDS(file=file.path(infodir, "panc.rna.rds"))

# load ATAC Seurat object
panc.atac <- readRDS(file=file.path(infodir, "panc.atac.rds"))

# perform cell clustering based on RNA data
my.dims <- 1:50
panc.rna %<>% FindNeighbors(reduction="mnn", dims=my.dims)
panc.rna %<>% FindClusters(resolution=0.7, verbose=T)

# Run UMAP dimensionality reduction
panc.rna %<>% RunUMAP(dims=my.dims, reduction="mnn", n.components=2, seed.use=42, 
                      n.neighbors=30, n.epochs=500)

# Set cell identity
panc.rna %<>% SetIdent(value="RNA_snn_res.0.7")

# Merge clusters
merge.clust <- c(0,4,0,2,0,1,3,1,3,2,1,2,2)
names(merge.clust) <- 0:12

final.clust <- as.vector(merge.clust[as.vector(Idents(panc.rna))])
names(final.clust) <- names(Idents(panc.rna))
final.clust <- factor(final.clust, levels=0:4)

# Add final clusters to meta data
panc.rna[["merged.cluster"]] <- final.clust

# Set final clusters
Idents(panc.rna) <- "merged.cluster"

# save Seurat object
saveRDS(panc.rna, file=file.path(infodir, "panc.rna.rds"))

# Run UMAP dimensionality reduction on ATAC data using the integrated embeddings
panc.atac %<>% RunUMAP(reduction = "integrated_lsi", dims = 2:50, 
                       n.components=2, seed.use=42, n.neighbors=30, n.epochs=500)

# save Seurat object
saveRDS(panc.atac, file=file.path(infodir, "panc.atac.rds"))

# Merge RNA and ATAC data into a single Seurat object
panc <- panc.rna
panc[["ATAC"]] <- panc.atac[['ATAC']]
# Add integrated_lsi dimensional reductions from integrated ATAC object
panc[["integrated_lsi"]] <- CreateDimReducObject(embeddings=Embeddings(panc.atac, reduction="integrated_lsi"),
                                                 key="IntegratedLSI_", assay="ATAC")
# Add UMAP dimensional reductions from integrated ATAC object
panc[["atac_integrated_umap"]] <- CreateDimReducObject(embeddings=Embeddings(panc.atac, reduction="umap"),
                                                       key="ATACIntegratedUMAP_", assay="ATAC")
# Rename UMAP dimensional reductions from RNA object
panc[["rna_umap"]] <- CreateDimReducObject(embeddings=Embeddings(panc, reduction="umap"),
                                           key="RNAUMAP_", assay="RNA")

# Compute a joint neighbor graph that represent both the gene expression and DNA accessibility measurements.
panc <- FindMultiModalNeighbors(
  object = panc,
  reduction.list = list("mnn", "integrated_lsi"), 
  dims.list = list(1:50, 2:50),
  verbose = TRUE
)

# Run joint UMAP dimensionality reduction
panc <- RunUMAP(
  object = panc,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
  seed.use=42,
  verbose = TRUE,
  n.neighbors=30, 
  n.epochs=500
)

# save Seurat object
saveRDS(panc, file.path(infodir, 'panc.rds'))
