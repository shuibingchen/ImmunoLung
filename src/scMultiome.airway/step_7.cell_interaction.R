# step_7.cell_interaction.R
# --- process multiome (RNA+ATAC) data from hPSC-derived immune-airway organoids ---
# step 7: Conduct cell-cell interactions analysis based on RNA data 
# Author: Tuo Zhang
# Date: 12/1/2024
# 

library(Seurat)
library(CellChat)
library(tidyverse)
library(patchwork)
library(BiocParallel)
library(magrittr)
library(future)
library(RColorBrewer)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)

# load functions
setwd(workdir)
source("my_functions.R")

# set parallelization in Seurat
plan("multicore", workers=6)
options(future.globals.maxSize=10*1024^3)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# load Seurat object
panc.rna <- readRDS(file.path(infodir, 'panc.rna.rds'))

# prepare a metadata sheet
rename.clusters <- c('0'='Ciliated', '1'='Basal-1', '2'='Fibroblast', '3'='Macrophage', '4'='Basal-2')

meta <- FetchData(panc.rna, vars=c('ident')) %>% mutate(ident=as.vector(ident)) %>% 
  mutate(ident=rename.clusters[ident]) %>% 
  dplyr::rename('labels'='ident')

# preprocess UMAP of hPSC-derived immune-airway organoids containing unstimulated (ARUM)
cellchat.arum <- my.preprocessing.cellchat(seurat.obj=panc.rna, sample.use=c('ARUM'), 
                                           clusters.use=0:4, metadata=meta, group.by='labels',
                                           CellChatDB.use=CellChatDB.use, 
                                           population.size=TRUE, infodir=infodir, 
                                           suffix='arum', 
                                           sources.use = 1:5, targets.use = 1:5, workers=6)

# preprocess UMAP of hPSC-derived immune-airway organoids containing proinflammatory macrophages (ARPM)
cellchat.arpm <- my.preprocessing.cellchat(seurat.obj=panc.rna, sample.use=c('ARPM'), 
                                           clusters.use=0:4, metadata=meta, group.by='labels',
                                           CellChatDB.use=CellChatDB.use, 
                                           population.size=TRUE, infodir=infodir, 
                                           suffix='arpm', 
                                           sources.use = 1:5, targets.use = 1:5, workers=6)

# merge cellchat objects
object.list <- list(ARUM = cellchat.arum, ARPM = cellchat.arpm)
cellchat.arum_arpm <- mergeCellChat(object.list, add.names = names(object.list))

# save object
saveRDS(cellchat.arum_arpm, file.path(infodir, 'cellchat.arum_arpm.rds'))
