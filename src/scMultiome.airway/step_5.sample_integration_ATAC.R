# step_5.sample_integration_ATAC.R
# --- process multiome (RNA+ATAC) data from hPSC-derived immune-airway organoids ---
# step 4: Integrate samples based on ATAC data
# Author: Tuo Zhang
# Date: 12/1/2024
# 

library(Seurat)
library(Signac)
library(tidyverse)
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

# load functions
setwd(workdir)
source("my_functions.R")

# load ATAC Seurat object list
seurat.obj.list <- readRDS(file.path(infodir, 'seurat.obj.list.rds'))

# load RNA Seurat object
panc.rna <- readRDS(file=file.path(infodir, "panc.rna.rds"))

# subset ATAC data, selecting cells that pass both RNA and ATAC QC filtering
for (sid in names(seurat.obj.list)){
  filtered.cells <- intersect(colnames(seurat.obj.list[[sid]]), colnames(panc.rna))
  seurat.obj.list[[sid]] <- subset(seurat.obj.list[[sid]], cells=filtered.cells)
}

# compute LSI
seurat.obj.list <- lapply(seurat.obj.list, FUN=my.process.ATAC, min.cutoff=20)

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = seurat.obj.list,
  anchor.features = rownames(seurat.obj.list[['ARUM']]),
  reduction = "rlsi",
  dims = 2:50
)

# integrate LSI embeddings
panc.atac <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = panc.atac[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
)

# save Seurat object
saveRDS(panc.atac, file=file.path(infodir, "panc.atac.rds"))
