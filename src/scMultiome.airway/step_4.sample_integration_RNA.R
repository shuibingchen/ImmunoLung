# step_4.sample_integration_RNA.R
# --- process multiome (RNA+ATAC) data from hPSC-derived immune-airway organoids ---
# step 4: Integrate samples based on RNA data
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

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# Set a random seed
set.seed(98)

# sample info
sample.info <- data.frame(SeqName=c("ARUM","ARPM"), 
                          Name=c("ARUM","ARPM"))
rownames(sample.info) <- sample.info$Name

# load Seurat object - RNA
panc.initial <- readRDS(file.path(infodir, 'panc.initial.rds'))

# load Seurat object - ATAC
panc.atac <- readRDS(file.path(infodir, 'panc.atac.rds'))

# fetch cells that pass ATAC QC filtering
atac.cells <- colnames(panc.atac)

# load (RNA-based) doublets detection result
df.results <- read.table(file.path(infodir, "doublet.results.txt"), header=T, check.names=F, stringsAsFactors=F, sep='\t')

# add doublets info to seurat object
panc.initial %<>% AddMetaData(metadata=df.results %>% dplyr::select('cellID','df.isdoublet') %>% column_to_rownames('cellID'))

# take intersecting cells that pass both RNA and ATAC QC filtering
singlet.cells <- rownames(subset(FetchData(panc.initial, vars=c('df.isdoublet')), df.isdoublet == 'Singlet'))
passing.filter.cells <- intersect(intersect(colnames(panc.initial), atac.cells), singlet.cells)

# create a new Seurat object containing both filtered RNA and ATAC data
panc.rna <- subset(panc.initial, cells=passing.filter.cells)

# create SingleCellExperiment object
sce.list <- list()
for (sample in rownames(sample.info)){
  print(sample)
  sce.list[[sample]] <- SingleCellExperiment(list(counts=panc.rna[["RNA"]]@counts[,rownames(subset(panc.rna@meta.data, orig.ident==sample))]))
}

# run a pre-clustering to avoid pooling together very different cells
# normalization will be performed for cells within each cluster
preclust.list <- lapply(sce.list, function(x) quickCluster(x=x, assay.type="counts"))

# normalize data by deconvolving size factors from cell pools
sce.list <- mapply(FUN=function(x,y) {computeSumFactors(x=x, min.mean=0.1, cluster=y)}, x=sce.list, y=preclust.list)

# compute normalized log-expression values
sce.list %<>% lapply(FUN=function(x) {logNormCounts(x)})

# rescale among donors (still need to update the normalized expressions in Seurat object)
rescaled.sce.list <- do.call(multiBatchNorm, sce.list)

# compute batch-corrected values across samples.
quick.corrected <- do.call(quickCorrect, c(rescaled.sce.list, list(PARAM=FastMnnParam(BSPARAM=BiocSingular::RandomParam()))))

# set column names
colnames(reducedDim(quick.corrected$corrected,"corrected")) = paste0("MNN_", 1:ncol(reducedDim(quick.corrected$corrected,"corrected")))

# Add normalized expression and MNN corrected data to Seurat object
panc.rna[["RNA"]]@data <- as(do.call(cbind, lapply(rescaled.sce.list, function(x) logcounts(x))) * log(2), "dgCMatrix")

# assign variable genes used for MNN correction
VariableFeatures(panc.rna) <- quick.corrected$hvgs

# add MNN correction results to Seurat object
panc.rna[["mnn"]] <- CreateDimReducObject(embeddings=reducedDim(quick.corrected$corrected, "corrected")[rownames(panc.rna@meta.data),], 
                                          key="MNN_", assay=DefaultAssay(panc.rna))

# save Seurat object
saveRDS(panc.rna, file=file.path(infodir, "panc.rna.rds"))
