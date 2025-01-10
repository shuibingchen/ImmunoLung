# step_7.cell_interaction.R
# --- process multiome (RNA+ATAC) data from hPSC-derived immune-alveolar organoids ---
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
rename.clusters <- c('0'='AT2', '1'='AT1', '2'='Basal-1', '3'='Macrophage', 
                     '4'='Basal-2', '5'='Basal-3', '6'='Fibroblast', '7'='Unknown')

meta <- FetchData(panc.rna, vars=c('ident')) %>% mutate(ident=as.vector(ident)) %>% 
  mutate(ident=rename.clusters[ident]) %>% 
  dplyr::rename('labels'='ident')

# preprocess UMAP of hPSC-derived immune-alveolar organoids containing unstimulated (AVUM)
cellchat.avum <- my.preprocessing.cellchat(seurat.obj=panc.rna, sample.use=c('AVUM'), 
                                           clusters.use=0:7, metadata=meta, group.by='labels',
                                           CellChatDB.use=CellChatDB.use, 
                                           population.size=TRUE, infodir=infodir, 
                                           suffix='avum', 
                                           sources.use = 1:8, targets.use = 1:8, workers=12)

# preprocess UMAP of hPSC-derived immune-alveolar organoids containing proinflammatory macrophages (AVPM)
cellchat.avpm <- my.preprocessing.cellchat(seurat.obj=panc.rna, sample.use=c('AVPM'), 
                                           clusters.use=0:7, metadata=meta, group.by='labels',
                                           CellChatDB.use=CellChatDB.use, 
                                           population.size=TRUE, infodir=infodir, 
                                           suffix='avpm', 
                                           sources.use = 1:8, targets.use = 1:8, workers=12)

# merge cellchat objects
object.list <- list(AVUM = cellchat.avum, AVPM = cellchat.avpm)
cellchat.avum_avpm <- mergeCellChat(object.list, add.names = names(object.list))

# save object
saveRDS(cellchat.avum_avpm, file.path(infodir, 'cellchat.avum_avpm.rds'))
