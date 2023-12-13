setwd("~/Projects/scCUT_TAG/")
library(GenomicRanges)
library(Seurat)
library(Signac)
library(ggplot2)
library(future)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(stringr)

counts <- Read10X_h5(filename = "PBMC/K27ac_PBMC_l1/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(file="PBMC/K27ac_PBMC_l1/outs/singlecell.csv",
                     header=TRUE,
                     row.names = 1)
chrom_assay <- CreateChromatinAssay(counts=counts,
                                    sep=c(":","-"),
                                    fragments = "PBMC/K27ac_PBMC_l1/outs/fragments.tsv.gz")
pbmc <- CreateSeuratObject(counts = chrom_assay,
                           assay="ATAC",
                           meta.data = metadata)


# extract gene annotations from EnsDb
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
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 1:50)

pbmc.rna <- readRDS("10x_scRNA/pbmc_10k_v3.rds")
pbmc.rna$tech <- "rna"

p1 <- DimPlot(pbmc, reduction = "umap") + NoLegend() + ggtitle("scCUT&TAG-seq")
p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))

transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc, 
                                        features = VariableFeatures(object = pbmc.rna), 
                                        reference.assay = "RNA", 
                                        query.assay = "ACTIVITY", 
                                        reduction = "cca")

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc.rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
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
  repel = TRUE) + NoLegend() + ggtitle('scCUT&TAG-seq')

plot1 + plot2
# output <- list(pbmc_rna= pbmc.rna,
#                pbmc_atac=pbmc)

predicted_id <- as.data.frame(pbmc$predicted.id)
predicted_id$barcode <- row.names(predicted_id)
colnames(predicted_id)[1] <- "predicted_id"

write.table(predicted_id,"atac_cuttag_comparison/cluster_compare/pbmc_cuttag_predicted_id.txt",sep="\t",
            row.names = F, quote = F)
