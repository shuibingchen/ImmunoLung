# step_3.sample_integration.R
# --- process scRNAseq data from hPSC-derived immune-alveolar organoids ---
# step 3: Integrate samples using MNN based corrections
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

# pattern for defining mitochondrial/ribosomal/virus genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"
cov2.pattern <- "^CoV2"

# Set a random seed
set.seed(98)

# sample info
sample.info <- data.frame(SeqName=paste('sc', 1:3, sep='_'),
                          Condition=c('Macrophage','Macrophage','293T'),
                          Infection=c('Mock','Cov2','Cov2'),
                          Name=c('AVMM','AVMS','AVTS'))
rownames(sample.info) <- sample.info$Name

# Load Seurat object
panc.initial <- readRDS(file.path(infodir, 'panc.initial.rds'))

# Load doublets result
df.results <- read.table(file.path(infodir, "doublet.results.txt"), 
                         header=T, check.names=F, stringsAsFactors=F, sep='\t')

# add doublets info to seurat object
panc.initial %<>% AddMetaData(metadata=df.results %>% 
                                dplyr::select('cellID','df.isdoublet') %>% 
                                column_to_rownames('cellID'))

# Remove doublets and create a new Seurat object for downstream analysis
singlet.cells <- rownames(subset(FetchData(panc.initial, vars=c('df.isdoublet')), 
                                 df.isdoublet == 'Singlet'))

# create a new Seurat object
panc <- subset(panc.initial, cells=singlet.cells)

# create SingleCellExperiment object
selected.samples <- rownames(sample.info)
sce.list <- list()
for (sample in selected.samples){
  print(sample)
  sce.list[[sample]] <- SingleCellExperiment(list(counts=panc[["RNA"]]@counts[,rownames(subset(panc@meta.data, orig.ident==sample))]))
}

# run a pre-clustering to avoid pooling together very different cells
# normalization will be performed for cells within each cluster
preclust.list <- lapply(sce.list, function(x) quickCluster(x=x, assay.type="counts"))

# normalize data by deconvolving size factors from cell pools
sce.list <- mapply(FUN=function(x,y) { computeSumFactors(x=x, min.mean=0.1, cluster=y) }, 
                   x=sce.list, y=preclust.list)

# compute normalized log-expression values
sce.list %<>% lapply(FUN=function(x) { logNormCounts(x) })

# rescale among donors
rescaled.sce.list <- do.call(multiBatchNorm, sce.list)

# discard mitochondrial, ribosomal and cov2 genes from integration
usable.genes <- grep(paste(mito.pattern, ribo.pattern, cov2.pattern, sep='|'), 
                     rownames(rescaled.sce.list[['AVMM']]), invert=T, value=T)

# compute batch-corrected values across samples.
quick.corrected <- do.call(quickCorrect, c(lapply(rescaled.sce.list, FUN=function(x) { x[usable.genes,] }), 
                                           list(PARAM=FastMnnParam(BSPARAM=BiocSingular::RandomParam()))))

# set column names
colnames(reducedDim(quick.corrected$corrected,"corrected")) <- paste0("MNN_", 1:ncol(reducedDim(quick.corrected$corrected,"corrected")))

# Add normalized expression and MNN corrected data to Seurat object
panc[["RNA"]]@data <- as(do.call(cbind, lapply(rescaled.sce.list, function(x) logcounts(x))) * log(2), "dgCMatrix")

# assign variable genes used for MNN correction
VariableFeatures(panc) <- quick.corrected$hvgs

# add MNN correction results to Seurat object
panc[["mnn"]] <- CreateDimReducObject(embeddings=reducedDim(quick.corrected$corrected, "corrected")[rownames(panc@meta.data),], 
                                      key="MNN_", assay=DefaultAssay(panc))

# save Seurat object
saveRDS(panc, file=file.path(infodir, "panc.rds"))
