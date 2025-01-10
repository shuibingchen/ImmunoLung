# step_1.qc.RNA.R
# --- process multiome (RNA+ATAC) data from hPSC-derived immune-airway organoids ---
# step 1: Load scRNA UMI counts data and perform cell QC
# Author: Tuo Zhang
# Date: 12/1/2024
# 

library(Seurat)
library(scater)
library(scran)
library(batchelor)
library(tidyverse)
library(scuttle)
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

# project name
project <- "yuling"

# pattern for defining mitochondrial/ribosomal genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"

# random seed
rseed <- 98
set.seed(rseed)

# load functions
setwd(workdir)
source("my_functions.R")

# sample info
sample.info <- data.frame(SeqName=c("ARUM","ARPM"), 
                          Name=c("ARUM","ARPM"))
rownames(sample.info) <- sample.info$Name

# load corrected UMI counts table per patient
raw.counts.list <- list()
for (k in 1:nrow(sample.info)){
  pid <- rownames(sample.info)[k]
  sid <- sample.info$SeqName[k]
  raw.counts.list[[k]] <- my.Read10X.GEX(file.path(sourcedir, sid, 'filtered_feature_bc_matrix'), pid)
}
names(raw.counts.list) <- rownames(sample.info)

# merge raw UMI counts tables
raw.counts.all <- my.MergeMatrix.v2(raw.counts.list)

# Initialize Seurat object
panc.initial <- CreateSeuratObject(counts=raw.counts.all, project=project, assay="RNA", 
                                   min.cells=0, min.features=0, 
                                   names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal genes per cell
panc.initial[["percent.mito"]] <- PercentageFeatureSet(panc.initial, pattern=mito.pattern)
panc.initial[["percent.ribo"]] <- PercentageFeatureSet(panc.initial, pattern=ribo.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc.initial@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc.initial@meta.data[,"orig.ident"])])
}
panc.initial %<>% AddMetaData(metadata=tmeta)

# calculate the number of cells pass filtering
qc.stats <- FetchData(panc.initial, vars=c('Name','nFeature_RNA','nCount_RNA','percent.mito'))
qc.stats <- qc.stats %>% 
  group_by(Name) %>% summarize_at('nFeature_RNA', list(nPre=length)) %>% 
  left_join(qc.stats %>% filter(nFeature_RNA > 300 & nFeature_RNA <= 7000 & 
                                  nCount_RNA > 600 & nCount_RNA <= 30000 & percent.mito < 20) %>% 
              group_by(Name) %>% summarize_at('nFeature_RNA', list(nPost=length)), by='Name') %>% 
  mutate(pRemain=round(nPost/nPre*100,1))
write.table(qc.stats, file.path(infodir,'qc.cell.number.txt'), quote=F, sep='\t', row.names=F)

# perform cell filtering
panc.initial %<>% subset(subset=nFeature_RNA > 300 & nFeature_RNA <= 7000 & 
                           nCount_RNA > 600 & nCount_RNA <= 30000 & percent.mito < 20)

# save Seurat object (initial)
saveRDS(panc.initial, file=file.path(infodir, "panc.initial.rds"))
