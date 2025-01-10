# step_5.differential_expression.R
# --- process scRNAseq data from hPSC-derived immune-alveolar organoids ---
# step 5: perform differential expression analysis
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
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# set a random seed
set.seed(98)

# Load Seurat object
panc <- readRDS(file.path(infodir, 'pand.rds'))

# define comparison groups
tinfo <- FetchData(panc, vars=c('ident','Name')) %>% 
  mutate(compare=paste0('C', ident, '-', Name))

# add to the seurat object
panc %<>% AddMetaData(metadata=tinfo[,'compare',drop=F])

# set compare clusters
Idents(panc) <- 'compare'

# all possible groups
all.groups <- c()
for (clust in 0:7){
  all.groups <- c(all.groups, paste0('C',clust,'-',c('AVMM','AVMS','AVTS')))
}

# get number of cells per group
nCells <- FetchData(panc, vars=c('ident')) %>% 
  dplyr::count(ident) %>% tidyr::complete(ident=all.groups, fill=list(n=0))

# perform DE analysis on AT2 cell cluster
for (clust in 0){
  # Immune-alveolar organoids exposed to SARS-CoV-2 vs. 
  # alveolar organoids co-cultured with 293T cells exposed to SARS-CoV-2
  my.DE.pair(panc, paste0('C',clust,'-AVMS'), paste0('C',clust,'-AVTS'), ntop=20, nCells=nCells,
             ttitle=paste0('C',clust,' AVM+S vs. AVT+S'), 
             no.mito=FALSE, mincells=10, min.pct=0, logfc.threshold=0, test.use='wilcox',
             cut.padj=0.1, cut.avglogfc=0, cutFC=0.6, cutP=20, 
             tfigdir=figdir, tinfodir=infodir)
}

# set original clusters back
Idents(panc) <- 'merged.cluster'

# load pathway gene sets
pathways <- read.delim(file.path(sourcedir, "pathways.txt"), header=T, sep='\t', 
                       check.names=T, stringsAsFactors=F)

# conduct gene set enrichment analysis
for (cluster in c(0)){
  # Immune-alveolar organoids exposed to SARS-CoV-2 vs. 
  # alveolar organoids co-cultured with 293T cells exposed to SARS-CoV-2
  cmp <- paste0('C',cluster,'-AVMS.vs.C',cluster,'-AVTS')
  # DE result file
  de.file <- file.path(infodir, paste('DE',cmp,'mincells_10.wilcox.min_pct_0.logfc_0.txt',sep='.'))
  # run GSEA
  my.run.GSEA(de.file=de.file, outdir=infodir, 
              suffix=paste('DE',cmp,'mincells_10.wilcox',sep='.'), 
              pathways=pathways, 
              minGSSize=10, maxGSSize=500)
}
