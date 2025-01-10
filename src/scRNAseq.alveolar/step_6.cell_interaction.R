# step_6.cell_interaction.R
# --- process scRNAseq data from hPSC-derived immune-alveolar organoids ---
# step 6: Conduct cell-cell interactions analysis between macrophages and AT2 cells 
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
panc <- readRDS(file=file.path(infodir, 'panc.rds'))

# remove cluster 7 (293T) since AVM+S does not have any 293 cells
panc.clean <- subset(panc, idents=c(0:6))

# prepare a metadata sheet
rename.clusters <- c('0'='AT2', '1'='Basal-1', '2'='AT1', '3'='Macrophage', '4'='Basal-2', 
                     '5'='Basal-3', '6'='Fibroblast')

meta <- FetchData(panc.clean, vars=c('ident')) %>% mutate(ident=as.vector(ident)) %>% 
  mutate(ident=rename.clusters[ident]) %>% 
  dplyr::rename('labels'='ident')

# preprocess immune-alveolar organoids exposed to mock (AVM+M)
cellchat.avmm <- my.preprocessing.cellchat(seurat.obj=panc.clean, sample.use=c('AVMM'), 
                                           clusters.use=0:6, metadata=meta, group.by='labels',
                                           CellChatDB.use=CellChatDB.use, 
                                           population.size=TRUE, infodir=infodir, 
                                           suffix='avmm', 
                                           sources.use = 1:7, targets.use = 1:7, workers=8)

# preprocess immune-alveolar organoids exposed to SARS-CoV-2 (MOI=0.1, AVM+S)
cellchat.avms <- my.preprocessing.cellchat(seurat.obj=panc.clean, sample.use=c('AVMS'), 
                                           clusters.use=0:6, metadata=meta, group.by='labels',
                                           CellChatDB.use=CellChatDB.use, 
                                           population.size=TRUE, infodir=infodir, 
                                           suffix='avms', 
                                           sources.use = 1:7, targets.use = 1:7, workers=8)

# merge cellchat objects
object.list <- list(AVMM = cellchat.avmm, AVMS = cellchat.avms)
cellchat.avmm_avms <- mergeCellChat(object.list, add.names = names(object.list))

# save object
saveRDS(cellchat.avmm_avms, file.path(infodir, 'cellchat.avmm_avms.rds'))
