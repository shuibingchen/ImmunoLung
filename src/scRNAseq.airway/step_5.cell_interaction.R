# step_5.cell_interaction.R
# --- process scRNAseq data from hPSC-derived immune-airway organoids ---
# step 5: Conduct cell-cell interactions analysis between macrophages and ciliated cells 
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

# remove cluster 6 (293T) since ARM+M/ARM+S have very few 293T cells
panc.clean <- subset(panc, idents=c(0:5))

# prepare a metadata sheet
rename.clusters <- c('0'='Basal-1', '1'='Ciliated', '2'='Macrophage', '3'='Basal-2', 
                     '4'='Virus-infected', '5'='Fibroblast')

meta <- FetchData(panc.clean, vars=c('ident')) %>% mutate(ident=as.vector(ident)) %>% 
  mutate(ident=rename.clusters[ident]) %>% 
  dplyr::rename('labels'='ident')

# preprocess immune-airway organoids exposed to mock (ARM+M)
cellchat.armm <- my.preprocessing.cellchat(seurat.obj=panc.clean, sample.use=c('ARMM'), 
                                           clusters.use=0:5, metadata=meta, group.by='labels',
                                           CellChatDB.use=CellChatDB.use, 
                                           population.size=TRUE, infodir=infodir, 
                                           suffix='armm', 
                                           sources.use = 1:6, targets.use = 1:6, workers=8)

# preprocess immune-airway organoids exposed to SARS-CoV-2 (MOI=0.1, AVM+S)
cellchat.arms <- my.preprocessing.cellchat(seurat.obj=panc.clean, sample.use=c('ARMS'), 
                                           clusters.use=0:5, metadata=meta, group.by='labels',
                                           CellChatDB.use=CellChatDB.use, 
                                           population.size=TRUE, infodir=infodir, 
                                           suffix='arms', 
                                           sources.use = 1:6, targets.use = 1:6, workers=8)

# merge cellchat objects
object.list <- list(ARMM = cellchat.armm, ARMS = cellchat.arms)
cellchat.armm_arms <- mergeCellChat(object.list, add.names = names(object.list))

# save object
saveRDS(cellchat.armm_arms, file.path(infodir, 'cellchat.armm_arms.rds'))
