setwd("~/Projects/scCUT_TAG/")
library(GenomicRanges)
library(Seurat)
library(Signac)
library(ggplot2)
library(future)
library(EnsDb.Hsapiens.v75)


peaks <- Read10X_h5("atac_pbmc_5k/atac_pbmc_5k_nextgem_filtered_peak_bc_matrix.h5")
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")

# create a gene activity matrix from the peak matrix and GTF, using chromosomes 1:22, X, and Y.
# Peaks that fall within gene bodies, or 2kb upstream of a gene, are considered
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(pbmc.atac) <- annotations

gene.activities <- GeneActivity(pbmc)
pbmc.atac[["ACTIVITY"]] <- pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)


