# my_functions.R
# customized functions for processing data and plotting
# Author: Tuo Zhang
# Date: 12/1/2024
# 

library(Seurat, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(Matrix, quietly=T, warn.conflicts=F)
library(CellChat, quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(pheatmap, quietly=T, warn.conflicts=F)
library(magrittr, quietly=T, warn.conflicts=F)
library(cowplot, quietly=T, warn.conflicts=F)
library(R.utils, quietly=T, warn.conflicts=F)
library(SingleCellExperiment, quietly=T, warn.conflicts=F)
library(tidyr, quietly=T, warn.conflicts=F)
library(tibble, quietly=T, warn.conflicts=F)
library(ggrepel, quietly=T, warn.conflicts=F)
library(ggplotify, quietly=T, warn.conflicts=F)
library(BiocParallel)
library(magrittr)
library(future)

# read in raw counts data from a given sample
my.Read10X <- function(tcountdir, tprefix=NULL){
  # directory exists?
  if (! dir.exists(tcountdir)){
    print(paste("input raw counts folder does NOT exist:", tcountdir, sep=" "))
    return(NULL)
  }
  # file exists?
  tmat.file <- paste(tcountdir, "matrix.mtx.gz", sep="/")
  tfnames.file <- paste(tcountdir, "features.tsv.gz", sep="/")
  tbnames.file <- paste(tcountdir, "barcodes.tsv.gz", sep="/")
  for (tf in c(tmat.file, tfnames.file, tbnames.file)){
    if (! file.exists(tf)){
      print(paste("input file does NOT exist:", tf))
      return(NULL)
    }
  }
  # extract counts matrix
  cat(paste("Loading UMI counts table from", tcountdir, "..."))
  tmat <- readMM(paste(tcountdir, "matrix.mtx.gz", sep="/"))
  tfnames <- read.delim(paste(tcountdir, "features.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  tbnames <- read.delim(paste(tcountdir, "barcodes.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  # update column names (cell ids)
  if (is.null(tprefix)){
    colnames(tmat) = tbnames$V1
  }
  else{
    colnames(tmat) = paste(tprefix, tbnames$V1, sep="_")
  }
  #rownames(tmat) = tfnames$V1
  # replace rowname (Ensembl id) by gene symbol
  # in case gene symbol is not unique, append the _EnsemblID after it
  # missing gene symbol will be replaced by EnsemblID
  #tsymbols <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(tmat), column="GENENAME", keytype="GENEID")
  ##tsymbols <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(tmat), column="SYMBOL", keytype="GENEID")
  rownames(tmat) <- uniquifyFeatureNames(ID=tfnames$V1, names=tfnames$V2)
  cat(" done.","\n")
  return(tmat)
}

# merge raw read counts table collected from multiple samples, in a more efficient way
my.MergeMatrix.v2 <- function(tmats){
  cat("Merge raw UMI counts ")
  tfunc <- function(x,y){
    tres <- cbind(x,y[rownames(x),])
    cat(".")
    return(tres)
  }
  tmerged <- Reduce(f=tfunc, x=tmats)
  # fill na with 0
  tmerged[is.na(tmerged)] <- 0
  cat(" done.")
  return(tmerged)
}

# function to compare pre-defined two groups of cells
# select top DE genes and make volcano plot
my.DE.pair <- function(tobj, cdtA, cdtB, ntop, nCells, ttitle, no.mito=TRUE, mincells=10,
                       min.pct=0.1, logfc.threshold=0.1, test.use='wilcox',
                       cut.padj=0.1, cut.avglogfc=0,
                       cutFC=0.6, cutP=20, tfigdir, tinfodir, tsvg=FALSE){
  # check number of cells in each condition
  nA <- as.numeric(nCells %>% filter(ident==cdtA) %>% dplyr::select(n))
  nB <- as.numeric(nCells %>% filter(ident==cdtB) %>% dplyr::select(n))
  print(paste0('Compare ',cdtA,' (',nA,' cells)',' to ',cdtB, ' (',nB,' cells).'))
  # enough cells?
  if (nA < mincells | nB < mincells){
    print('Too few cells to perform DE.')
    return(NULL)
  }
  # run DE
  de.results <- FindMarkers(tobj, ident.1=cdtA, ident.2=cdtB, min.pct=min.pct, logfc.threshold=logfc.threshold, test.use=test.use)
  # write to file
  out.file.prefix <- paste('DE',cdtA,'vs',cdtB,paste0('mincells_',mincells),
                           test.use,paste0('min_pct_',min.pct),paste0('logfc_',logfc.threshold), sep='.')
  de.file <- file.path(tinfodir, paste(out.file.prefix, 'txt', sep='.'))
  write.table(de.results, de.file, quote=F, sep='\t', col.names=NA)
  # select top 10 up/down genes to label in the plot
  up.genes <- rownames(subset(de.results, p_val_adj < cut.padj & avg_log2FC > cut.avglogfc))
  down.genes <- rownames(subset(de.results, p_val_adj < cut.padj & avg_log2FC < -cut.avglogfc))
  top.up.genes <- head(up.genes, ntop)
  top.down.genes <- head(down.genes, ntop)
  # ignore mitochondiral genes?
  if (no.mito){
    top.up.genes <- head(grep(mito.pattern, up.genes, value=T, invert=T), ntop)
    top.down.genes <- head(grep(mito.pattern, down.genes, value=T, invert=T), ntop)
  }
  # make volcano plot
  g <- myVolcanoPlot(de.file, tx='p_val_adj', ty='avg_log2FC', tcutFC=cutFC, tcutP=cutP, tlabel.genes.up=top.up.genes, tlabel.genes.down=top.down.genes,
                     txlabel='LogFoldChange', tylabel='-LogP-value', tupper=300, talpha=0.8,
                     tcolor.up='firebrick3', tcolor.down='steelblue3', tcolor.other='gray60')
  g <- g + ggtitle(ttitle) + theme(plot.title=element_text(color='black',face='bold',size=20,hjust=0.5))
  ggsave(file.path(tfigdir,paste('Volcano',out.file.prefix,'png',sep='.')), plot=g, width=8, height=7, dpi=300)
  if (tsvg){
    ggsave(file.path(tfigdir,paste('Volcano',out.file.prefix,'svg',sep='.')), plot=g, width=8, height=7)
  }
  return(de.results)
}

my.run.GSEA <- function(de.file, outdir, suffix, pathways=pathways, 
                        minGSSize=10, maxGSSize=500){
  # load DE results
  de.data <- read.table(de.file, header=T, check.names=F, stringsAsFactors=F, sep='\t', row.names=1)
  print(dim(de.data))
  
  # prepare input for GSEA (by Symbol)
  clean.de.data <- de.data %>% rownames_to_column('Symbol') %>% 
    arrange(desc(avg_log2FC))
  symbol.genes.gsea <- clean.de.data %>% pull(avg_log2FC)
  names(symbol.genes.gsea) <- clean.de.data %>% pull('Symbol')
  
  # run GSEA on pathway gene sets (output all pathways)
  gse <- GSEA(geneList=symbol.genes.gsea,
              TERM2GENE=pathways, 
              minGSSize=minGSSize, 
              maxGSSize=maxGSSize, 
              pvalueCutoff=1, 
              pAdjustMethod='BH', 
              seed = TRUE, 
              verbose=verbose)
  # save
  saveRDS(gse, file.path(outdir, paste('gse',suffix,'rds', sep='.')))
}

# make volcano plot on a set of DE genes
myVolcanoPlot <- function(tdefile, tx='p_val', ty='avg_logFC', tcutFC=0.25, tcutP=20, tlabel.genes.up=NULL, tlabel.genes.down=NULL, txlabel='LogFoldChange', tylabel='-LogP-value', tupper=300, talpha=0.8, tcolor.up='firebrick3', tcolor.down='steelblue3', tcolor.other='gray60'){
  # read in DE results
  tde <- read.table(tdefile, sep='\t', header=T, check.names=F, stringsAsFactors=F, row.names=1)
  # arrange data
  tdataToPlot <- tde[,c(tx,ty)]
  tdataToPlot$gene <- rownames(tde)
  colnames(tdataToPlot) <- c('pval','logFC','gene')
  # log transform p-value
  tdataToPlot$logPval <- -log10(tdataToPlot$pval)
  # label and color
  if (is.null(tlabel.genes.up) | is.null(tlabel.genes.down)){
    tdataToPlot$color <- with(tdataToPlot, ifelse(logPval > tcutP & logFC > tcutFC, 'up', ifelse(logPval > tcutP & logFC < -tcutFC, 'down', 'other')))
    tdataToPlot$label <- with(tdataToPlot, ifelse(color %in% c('up','down'), gene, ''))
  } else{
    tdataToPlot$color <- with(tdataToPlot, ifelse(gene %in% tlabel.genes.up, 'up', ifelse(gene %in% tlabel.genes.down, 'down', 'other')))
    tdataToPlot$label <- with(tdataToPlot, ifelse(gene %in% c(tlabel.genes.up, tlabel.genes.down), gene, ''))
  }
  # any 0 p-values thus Inf logPval? modify to the upperlimit
  tdataToPlot$logPval[tdataToPlot$logPval > tupper] <- tupper
  # plot
  tg <- ggplot(tdataToPlot, aes(x=logFC, y=logPval, color=color, label=label))
  tg <- tg + geom_point(shape=19, size=2, alpha=talpha)
  tg <- tg + scale_color_manual(values=c('up'=tcolor.up,'down'=tcolor.down, 'other'=tcolor.other))
  tg <- tg + geom_text_repel()
  #tg <- tg + geom_vline(xintercept=tcutFC, linetype='dotted', color='gray75') + geom_vline(xintercept=-tcutFC, linetype='dotted', color='gray75')
  #tg <- tg + geom_hline(yintercept=tcutP, linetype="dotted", color="gray75")
  tg <- tg + xlab(txlabel) + ylab(tylabel)
  tg <- tg + theme_classic()
  tg <- tg + theme(legend.position='none')
  tg <- tg + theme(axis.title=element_text(size=18, color='black'), axis.text=element_text(size=16, color='black'))
  return(tg)
}

# pre-processing cells in a particular group
my.preprocessing.cellchat <- function(seurat.obj, sample.use, clusters.use, metadata, group.by='labels',
                                      CellChatDB.use=CellChatDB.use, population.size=TRUE, infodir=infodir, 
                                      suffix='', sources.use = c(3,7), targets.use = c(3,7), workers=4){
  # define cells used for running CellChat
  tcell.use <- rownames(subset(FetchData(seurat.obj, vars=c('ident', 'Name')), 
                               Name %in% sample.use & ident %in% clusters.use))
  
  # prepare input data for running CellChat
  tdata.input <- seurat.obj[['RNA']]@data[, tcell.use]
  tmeta <- metadata[tcell.use, , drop=F]
  
  # create a CellChat object
  cellchat <- createCellChat(object=tdata.input, meta=tmeta, group.by=group.by)
  print(cellchat)
  
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  # Preprocessing the expression data for cell-cell communication analysis
  # subset the expression data of signaling genes for saving computation cost
  print('subsetData')
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

  print('identifyOverExpressedGenes')
  cellchat <- identifyOverExpressedGenes(cellchat)
  print('identifyOverExpressedInteractions')
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  # cellchat <- projectData(cellchat, PPI.human)
  
  # Compute the communication probability and infer cellular communication network
  print('computeCommunProb')
  cellchat <- computeCommunProb(cellchat, population.size=population.size)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Extract the inferred cellular communication network as a data frame
  df.net <- subsetCommunication(cellchat)
  write.table(df.net, file=file.path(infodir, paste('cell-cell.communication',suffix,'all.txt',sep='.')), quote=F, sep='\t', row.names=F)
  
  # interactions between beta and macrophage
  print(subsetCommunication(cellchat, sources.use = sources.use, targets.use = targets.use))
  
  # Infer the cell-cell communication at a signaling pathway level
  # CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.
  # NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
  print('computeCommunProbPathway')
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate the aggregated cell-cell communication network
  print('aggregateNet')
  cellchat <- aggregateNet(cellchat)
  
  # Compute the network centrality scores
  print('netAnalysis_computeCentrality')
  future::plan("sequential")
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  future::plan("multicore", workers = workers)
  
  # save CellChat object
  saveRDS(cellchat, file.path(infodir, paste('cellchat',suffix,'rds',sep='.')))
  
  return(cellchat)
}

# plot cells colored by an annotation (DimPlot)
myDimPlot <- function(tobj, treduct, tcate, torder=NULL, tsuffix, tcells=NULL, tcolor=NULL, 
                      tlabel=FALSE, tsplit=FALSE, txlim=NULL, tylim=NULL,
                      tncol=4, tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20){
  tdataToPlot <- data.frame()
  tlabel.pos <- data.frame()
  tg <- ggplot()
  if (treduct == "tsne"){
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("tSNE_1","tSNE_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      missed.categroy <- setdiff(unique(tdataToPlot$Category), torder)
      if (length(missed.categroy) > 0){
        tdataToPlot <- rbind(tdataToPlot, data.frame(UMAP_1=NA, UMAP_2=NA, Category=missed.categroy))
      }
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=tSNE_1, y=tSNE_2, color=Category))
    tg <- tg + ggtitle(paste("tSNE","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(tSNE_1, tSNE_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  } else if (treduct == "umap") {
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("UMAP_1","UMAP_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      missed.categroy <- setdiff(unique(tdataToPlot$Category), torder)
      if (length(missed.categroy) > 0){
        tdataToPlot <- rbind(tdataToPlot, data.frame(UMAP_1=NA, UMAP_2=NA, Category=missed.categroy))
      }
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=UMAP_1, y=UMAP_2, color=Category))
    tg <- tg + ggtitle(paste("UMAP","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(UMAP_1, UMAP_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  }
  if (! is.null(txlim)){
    tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha)
  if (! is.null(tcolor)){
    tg <- tg + scale_color_manual(values=tcolor)
  }
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (tsplit == TRUE){
    tg <- tg + facet_wrap(~Category, ncol=tncol)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Category), color="black")
  }
  return(tg)
}

myDimPlot3 <- function(tobj, treduct, tgroup_by, tgroup_order=NULL, thighlight=NULL, 
                       tsuffix, tcells=NULL, tcolor=NULL, tlabel=FALSE, tsplit_by=NULL, 
                       tsplit_order=NULL, txlim=NULL, tylim=NULL,
                       tncol=1, tptshape=19, tptsize=2, talpha=0.6, tltsize=18, 
                       tatlsize=20, tlbsize=2){
  # set coordinates variable name
  vars.reduct <- c("UMAP_1","UMAP_2")
  if (treduct == "tsne"){
    vars.reduct <- c("tSNE_1","tSNE_2")
  }
  # extract coordinates + group
  tdataToPlot <- FetchData(tobj, cells=tcells, vars=c(vars.reduct, tgroup_by))
  colnames(tdataToPlot) <- c("Dim_1","Dim_2","Group")
  # update group order if available
  if (!is.null(tgroup_order)){
    # add fake rows to make sure each group is considered
    tmiss.cates <- setdiff(tgroup_order, unique(tdataToPlot$Group))
    if (length(tmiss.cates) > 0){
      tdataToPlot <- rbind(tdataToPlot, data.frame(Dim_1=NA, Dim_2=NA, Group=tmiss.cates))
    }
    # reorder categories
    tdataToPlot$Group <- factor(tdataToPlot$Group, levels=tgroup_order)
  }
  # extract split
  if (! is.null(tsplit_by)){
    tsp <- FetchData(tobj, cells=tcells, vars=c(tsplit_by))
    tdataToPlot$Split <- tsp[rownames(tdataToPlot), c(tsplit_by)]
    # update split order if available
    if (! is.null(tsplit_order)){
      tdataToPlot$Split <- factor(tdataToPlot$Split, levels=tsplit_order)
    }
  }
  # reorder cells that needs to highlight (draw those cell points later)
  if (!is.null(thighlight)){
    tdataToPlot <- rbind(subset(tdataToPlot, ! Group %in% thighlight), subset(tdataToPlot, Group %in% thighlight))
  }
  # prepare group labeling
  tlabel.pos <- aggregate(cbind(Dim_1, Dim_2) ~ Group, data=tdataToPlot, FUN=median)
  colnames(tlabel.pos) <- c("Group","X","Y")
  # plot
  tg <- ggplot(tdataToPlot, aes(x=Dim_1, y=Dim_2, color=Group))
  tg <- tg + ggtitle(paste(toupper(treduct),"plots","by",tsuffix, sep=" "))
  # set range on coordinates
  if (! is.null(txlim)){
    tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha)
  if (! is.null(tcolor)){
    tg <- tg + scale_color_manual(values=tcolor)
  }
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (! is.null(tsplit_by)){
    tg <- tg + facet_wrap(~Split, ncol=tncol)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Group), color="black", size=tlbsize)
  }
  return(tg)
}

my.dot.plot <- function(seurat.obj, genes, cluster.id=0, name.order=c('AVMM','AVMC')){
  # extract subset cell ids
  tcells <- FetchData(seurat.obj, vars=c('ident', 'Name')) %>% 
    rownames_to_column('cellID') %>% 
    dplyr::filter(ident %in% !!cluster.id, Name %in% !!name.order) %>% 
    pull('cellID')
  # subset Seurat object
  tseu <- subset(seurat.obj, cells=tcells)
  # order samples
  tseu$Name <-  factor(tseu$Name, levels=name.order)
  # dot plot
  return(DotPlot(tseu, features=genes, group.by='Name') + coord_flip())
}
