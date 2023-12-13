setwd("~/Projects/scCUT_TAG/")
library(GenomicRanges)
library(Seurat)
library(Signac)
library(ggplot2)
library(future)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(stringr)

# input macs2 peaks
macs_peak <- read.table("PBMC/macs_analysis_l1/macs_out/PBMC_1_peaks.narrowPeak",sep="\t",header=F)

target_chr <- c("chr1","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18",
                "chr19","chr20","chr21","chr22","chrX")

# only keep ch1-chr22,chrX
macs_peak <- macs_peak[which(macs_peak$V1 %in% target_chr),]

macs_peak <- macs_peak[,c(1:3)]
colnames(macs_peak) <- c("chr","start","end")

# convert to genomic ranges
macs_gr <- makeGRangesFromDataFrame(macs_peak)

#------------------------------------
# quantify peaks 
#------------------------------------
#counts <- Read10X_h5(filename = "PBMC/K27ac_PBMC_l1/outs/filtered_peak_bc_matrix.h5")
# metadata <- read.csv(file="PBMC/K27ac_PBMC_l1/outs/singlecell.csv",
#                      header=TRUE,
#                      row.names = 1)

# create fragments objects
macs_frag <- CreateFragmentObject(path = "PBMC/K27ac_PBMC_l1/outs/fragments.tsv.gz")

# quantify peaks 
macs_count <- FeatureMatrix(
  fragments = macs_frag,
  features = macs_gr
)

# create Seurat object
pbmc_assay <- CreateChromatinAssay(counts=macs_count,
                                   sep=c(":","-"),
                                   fragments = macs_frag,
                                   min.cells = 10,
                                   min.features = 200)
pbmc <- CreateSeuratObject(counts = pbmc_assay,
                               assay="ATAC")

#--------------------------------------
# extract gene annotations from EnsDb
#--------------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# seqlevels(annotations)
# head(Fragments(pbmc)[[1]])

# change seqlevels from * to chr*
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), 
                           pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels

# add the gene information to the object
Annotation(pbmc) <- annotations

gene.activities <- GeneActivity(pbmc)

pbmc[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
pbmc$tech <- "atac"


DefaultAssay(pbmc) <- "ACTIVITY"
pbmc <- FindVariableFeatures(pbmc)
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)

DefaultAssay(pbmc) <- "ATAC"
VariableFeatures(pbmc) <- names(which(Matrix::rowSums(pbmc) > 100))
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
#pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 1:50)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 1:50)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 1:50)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()

#--------------------------------------
# get scRNA-seq
#--------------------------------------
pbmc.rna <- readRDS("10x_scRNA/pbmc_10k_v3.rds")
pbmc.rna$tech <- "rna"

p1 <- DimPlot(pbmc, reduction = "umap") + NoLegend() + ggtitle("scCUT&TAG-seq (Macs2)")
p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))

#--------------------------------------
# mapping scRNA and scATAC
#--------------------------------------
mem.maxVSize(vsize = Inf)
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc, 
                                        features = VariableFeatures(object = pbmc.rna), 
                                        reference.assay = "RNA", 
                                        query.assay = "ACTIVITY",
                                        reduction = "cca")


predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc.rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 1:50
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

plot1 <- DimPlot(
  object = pbmc.rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scCUT&TAG-seq (Peaks)')

plot1 + plot2
# output <- list(pbmc_rna= pbmc.rna,
#                pbmc_atac=pbmc)

predicted_id <- as.data.frame(pbmc$predicted.id)
predicted_id$barcode <- row.names(predicted_id)
colnames(predicted_id)[1] <- "predicted_id"

write.table(predicted_id,"atac_cuttag_comparison/cluster_compare/pbmc_cuttag_macs2_predicted_id.txt",sep="\t",
            row.names = F, quote = F)
