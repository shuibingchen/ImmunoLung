# ImmunoLung

This repository contains the codes necessary to perform the analysis in our 
manuscript [Human Immuno-Lung Organoid Model to Study Macrophage-Mediated Lung Cell Senescence Upon SARS-CoV-2 Infection](https://github.com/shuibingchen/ImmunoLung), 
as described in the methods and main text.

### Input data

The single cell data were generated with the 10X Genomics kits and
pre-processed using 10X cellranger pipeline. The data are available in the
GEO database with accession#[XXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=XXX)

### Requirements

The following R packages are needed

- Seurat
- Signac
- CellChat
- clusterProfiler
- scran
- scater
- batchelor
- tidyverse
- pheatmap
- RColorBrewer
- BiocParallel
- magrittr
- ggrepel
- enrichplot
- patchwork
- org.Hs.eg.db
- EnsDb.Hsapiens.v86

