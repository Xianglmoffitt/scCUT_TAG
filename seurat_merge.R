setwd("~/Projects/scCUT_TAG/")
library(GenomicRanges)
library(Seurat)
library(Signac)
library(ggplot2)
library(future)


target_chr <- c("chr1","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18",
                "chr19","chr20","chr21","chr22","chrX")

atac_peak <- read.table(
  file = "atac_pbmc_5k/atac_pbmc_5k_nextgem_peaks.bed",
  col.names = c("chr", "start", "end")
)

cuttag_peak <- read.table(
  file = "PBMC/K27ac_PBMC_l1/outs/peaks.bed",
  col.names = c("chr", "start", "end")
)

atac_peak <- atac_peak[which(atac_peak$chr %in% target_chr),]
cuttag_peak <- cuttag_peak[which(cuttag_peak$chr %in% target_chr),]

atac_peak_gr <- makeGRangesFromDataFrame(atac_peak)
cuttag_peak_gr <- makeGRangesFromDataFrame(cuttag_peak)

# combine peaks
combined_peak <- reduce(x=c(atac_peak_gr,cuttag_peak_gr))

# load metadata
atac_md <- read.table(
  file = "atac_pbmc_5k/atac_pbmc_5k_nextgem_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

cuttag_md <- read.table(
  file = "PBMC/K27ac_PBMC_l1/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]


# compute hash
frags_atac <- CreateFragmentObject(
  path = "atac_pbmc_5k/atac_pbmc_5k_nextgem_fragments.tsv.gz",
  cells = rownames(atac_md)
)

frags_cuttag <- CreateFragmentObject(
  path = "PBMC/K27ac_PBMC_l1/outs/fragments.tsv.gz",
  cells = rownames(cuttag_md)
)

# quantify peaks
atac_counts <- FeatureMatrix(
  fragments = frags_atac,
  features = combined_peak,
  cells = rownames(atac_md)
)

cuttag_counts <- FeatureMatrix(
  fragments = frags_cuttag,
  features = combined_peak,
  cells = rownames(cuttag_md)
)

# creat objects
atac_assay <- CreateChromatinAssay(atac_counts, fragments = frags_atac)
atac <- CreateSeuratObject(atac_assay, assay = "ATAC", meta.data=atac_md)

cuttag_assay <- CreateChromatinAssay(cuttag_counts, fragments = frags_cuttag)
cuttag <- CreateSeuratObject(cuttag_assay, assay = "ATAC", meta.data=cuttag_md)

# merge object
atac$dataset <- 'pbmc_atac'
cuttag$dataset <- 'pbmc_cuttag'

combined <- merge(
  x = atac,
  y = list(cuttag),
  add.cell.ids = c("atac", "cuttag")
)
combined[["ATAC"]]


combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)


