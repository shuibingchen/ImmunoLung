# step_2.qc.ATAC.R
# --- process multiome (RNA+ATAC) data from hPSC-derived immune-alveolar organoids ---
# step 2: Load scATAC UMI counts data and perform cell QC
# Author: Tuo Zhang
# Date: 12/1/2024
# 

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(patchwork)
library(future)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(GenomicRanges)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# project name
project <- "yuling"

# random seed
rseed <- 98
set.seed(rseed)

# load functions
setwd(workdir)
source("my_functions.R")

# set parallelization in Seurat
plan("multicore", workers=6)
options(future.globals.maxSize=10*1024^3)

# ATAC blacklist regions
blacklist.regions <- blacklist_hg38_unified

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to USCS style since the chromsomes start with 'chr'
seqlevels(annotations) <- ifelse(seqlevels(annotations) == 'MT', 'chrM', paste0('chr', seqlevels(annotations)))
genome(annotations) <- "hg38"

# sample info
sample.info <- data.frame(SeqName=c("AVUM","AVPM"), 
                          Name=c("AVUM","AVPM"))
rownames(sample.info) <- sample.info$Name

# load peak sets
gr.list <- list()

for (k in 1:nrow(sample.info)){
  pid <- sample.info$Name[k]
  sid <- sample.info$SeqName[k]
  print(sid)
  peaks <- read.table(file=file.path(sourcedir, sid, 'atac_peaks.bed'), 
                      col.names = c("chr", "start", "end"))
  gr.list[[pid]] <- makeGRangesFromDataFrame(peaks)
}

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- Signac::reduce(c(gr.list[[1]], gr.list[[2]]))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

# Create Seurat object for each sample
seurat.obj.list <- list()

for (k in 1:nrow(sample.info)){
  pid <- sample.info$Name[k]
  sid <- sample.info$SeqName[k]
  print(sid)
  # load metadata
  metadata <- read.table(file=file.path(sourcedir, sid, 'per_barcode_metrics.csv'), 
                         header=T, check.names=F, stringsAsFactors=F, sep=',', row.names=1)
  # select columns
  metadata <- metadata %>% dplyr::select(c('gex_barcode','atac_barcode','is_cell',
                                           'atac_fragments','atac_TSS_fragments',
                                           'atac_peak_region_fragments',
                                           'atac_peak_region_cutsites'))
  # fetch valid cells
  metadata <- metadata %>% dplyr::filter(is_cell == 1)
  cells <- rownames(metadata)
  # create fragment objects
  frags <- CreateFragmentObject(path=file.path(sourcedir, sid, 'atac_fragments.tsv.gz'), 
                                cells=cells)
  # quantify peaks
  # create a matrix of peaks x cell
  counts <- FeatureMatrix(fragments=frags, features=combined.peaks, cells=cells)
  # create object
  chrom.assay <- CreateChromatinAssay(counts=counts, fragments=frags)
  seurat.obj.list[[pid]] <- CreateSeuratObject(counts=chrom.assay, assay = "ATAC", meta.data=metadata)
  # add gene annotations
  Annotation(seurat.obj.list[[pid]]) <- annotations
}

# calculate cell QC metrics: 
# nucleosome signal score, TSS enrichment score, fraction of reads in peaks, blacklist ratio
seurat.obj.list <- lapply(seurat.obj.list, FUN=my.compute.QC.metrics)

# filter cells based on QC metrics
seurat.obj.list <- lapply(seurat.obj.list, FUN=my.filter.cells, 
                          nCount.min=700, nCount.max=25000, pct.rip.min=50, 
                          bl.frac.max=0.01, nuc.sig.max=1, tss.min=3)

# add information to identify dataset of origin
# add sample label to cell ids (AAACAGCCAATGCCCG-1 --> AVUM_AAACAGCCAATGCCCG-1)
for (k in 1:nrow(sample.info)){
  pid <- sample.info$Name[k]
  seurat.obj.list[[pid]]$dataset <- pid
  seurat.obj.list[[pid]] %<>% RenameCells(add.cell.id = pid)
}

# save seurat object list
saveRDS(seurat.obj.list, file.path(infodir, 'seurat.obj.list.rds'))
