setwd("~/Projects/scCUT_TAG/")
library(GenomicRanges)

# peak <- read.table("peaks.bed",sep="\t",header=F)
# peak$peak_name <- paste(peak$V1,peak$V2,peak$V3,sep="_")
# peak$strand <- rep(".",nrow(peak))
# write.table(peak,"../peak_tmp.gff",quote = F,col.names = F,row.names = F,sep="\t")


# cell_bar <- read.table("filtered_peak_bc_matrix/barcodes.tsv",sep="\t",header=F)
# 
# for (i in 1:nrow(cell_bar)){
#   out_name <- paste("cell",i,sep="_")
#   cell_bc <- cell_bar$V1[i]
#   write.table(cell_bc,paste0("../cell_bc/",out_name),quote=F,row.names = F,col.names = F,sep="\t")
# }

# extract example SE regions, chr1:153944937-153991211

peaks <- read.table("atac_pbmc_5k/filtered_peak_bc_matrix/peaks.bed",sep="\t",header=F)
bc_cell <- read.table("atac_pbmc_5k/filtered_peak_bc_matrix/barcodes.tsv",sep="\t",header=F)
count_matrix_atac_5k <- read.table("atac_pbmc_5k/filtered_peak_bc_matrix/matrix.mtx",sep=" ",header=F,
                           comment.char = "%")

# all cell
summary(count_matrix_atac_5k$V3[-1])
summary(count_matrix_cut_tag$V3[-1])

# one cell
summary(count_matrix_atac_5k$V3[which(count_matrix_atac_5k$V1==48742)])
summary(count_matrix_cut_tag$V3[which(count_matrix_atac_5k$V1==48742)])

shosen_peaks <- which(peaks$V1=="chr1"& peaks$V2>=153944937 &peaks$V3<=153991211)

sub_matrix <- count_matrix[which(count_matrix$V1 %in% shosen_peaks),]
sub_bc <- unique(sub_matrix$V2)

example_summary <- data.frame()
for (i in 1:length(sub_bc)) {
  example_m_tmp <- sub_matrix[which(sub_matrix$V2==sub_bc[i]),]
  summary_tmp <- data.frame(bc=bc_cell[sub_bc[i],],
                            peak_n=nrow(example_m_tmp),
                            umi_sum=sum(example_m_tmp$V3))
  example_summary <- rbind(example_summary,summary_tmp)
}

example_final <- example_summary[order(-example_summary$peak_n,-example_summary$umi_sum),]

write.table(example_final,"../se_example_summary_table.txt",sep="\t",row.names = F,quote=F)


#chr1:234599582_234620744

#-----------------------------------
# extract cluster barcorde
#-----------------------------------

cluster_df <- read.table("analysis/clustering/graphclust/clusters.csv",sep=",",header=T)
cluster_id <- unique(cluster_df$Cluster)
for (c in 1:length(cluster_id)) {
  tmp <- cluster_df$Barcode[which(cluster_df$Cluster==cluster_id[c])]
  outname <- paste0("../cell_clusters/cluster_",c,".txt")
  write.table(tmp,outname,quote=F,row.names = F,col.names = F)
}


#------------------------------------------------
# macs2 and cellranger peak comparison
#------------------------------------------------
target_chr <- c("chr1","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18",
                "chr19","chr20","chr21","chr22","chrX")
macs_peak <- read.table("PBMC/macs_analysis_l1/macs_out/PBMC_1_peaks.narrowPeak",sep="\t")
cellranger_peak <- read.table("PBMC/K27ac_PBMC_l1/outs/peaks.bed",sep="\t")
# only keep ch1-chr22,chrX
macs_peak <- macs_peak[which(macs_peak$V1 %in% target_chr),]
cellranger_peak <- cellranger_peak[which(cellranger_peak$V1 %in% target_chr),]

# add names
cellranger_peak$peak_name <- paste(cellranger_peak$V1,cellranger_peak$V2,
                                   cellranger_peak$V3,sep="_")
#atac_peak <- read.table("atac_pbmc_5k/atac_pbmc_5k_nextgem_peaks.bed",sep="\t")

# make granges
macs_peak_gr <- GRanges(
  seqnames=macs_peak$V1,
  ranges=IRanges(macs_peak$V2,macs_peak$V3,names=macs_peak$V4)
)

cellranger_peak_gr <- GRanges(
  seqnames=cellranger_peak$V1,
  ranges=IRanges(cellranger_peak$V2,cellranger_peak$V3,names=cellranger_peak$peak_name)
)

# overlap peaks
overlap_peaks <- findOverlaps(macs_peak_gr,cellranger_peak_gr)

# unique peaks
macs_peak_unique <- unique(macs_peak[-queryHits(overlap_peaks),])
cellranger_peak_unique <- unique(cellranger_peak[-subjectHits(overlap_peaks),])

overlap_df <- cbind(cellranger_peak[subjectHits(overlap_peaks),],macs_peak[queryHits(overlap_peaks),c(1:4)])

write.table(macs_peak_unique[,c(1:4)],"peak_compare/macs_peak_unique.txt",sep = "\t",
            quote = F,row.names = F,col.names = F)
write.table(cellranger_peak_unique,"peak_compare/cellranger_peak_unique.txt",sep = "\t",
            quote = F,row.names = F,col.names = F)
write.table(overlap_df,"peak_compare/overlap_df.txt",sep = "\t",
            quote = F,row.names = F,col.names = F)




#------------------------------------------------
# macs2 and cellranger SE comparison
#------------------------------------------------
target_chr <- c("chr1","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18",
                "chr19","chr20","chr21","chr22","chrX")
macs_se <- read.table("PBMC/macs_analysis_l1/rose_out/PBMC_1_peaks_Gateway_SuperEnhancers.bed",sep="\t")
cellranger_se <- read.table("PBMC/K27ac_PBMC_l1/rose_out/peak_l1_Gateway_SuperEnhancers.bed",sep="\t")
# only keep ch1-chr22,chrX
macs_se <- macs_se[which(macs_se$V1 %in% target_chr),]
cellranger_se <- cellranger_se[which(cellranger_se$V1 %in% target_chr),]

# add names
cellranger_se$peak_name <- paste(cellranger_se$V1,cellranger_se$V2,
                                   cellranger_se$V3,sep="_")
#atac_peak <- read.table("atac_pbmc_5k/atac_pbmc_5k_nextgem_peaks.bed",sep="\t")

# make granges
macs_se_gr <- GRanges(
  seqnames=macs_se$V1,
  ranges=IRanges(macs_se$V2,macs_se$V3,names=macs_se$V4)
)

cellranger_se_gr <- GRanges(
  seqnames=cellranger_se$V1,
  ranges=IRanges(cellranger_se$V2,cellranger_se$V3,names=cellranger_se$peak_name)
)

# overlap peaks
overlap_se <- findOverlaps(macs_se_gr,cellranger_se_gr)

# unique peaks
macs_peak_unique <- unique(macs_se[-queryHits(overlap_se),])
cellranger_peak_unique <- unique(cellranger_se[-subjectHits(overlap_se),])

overlap_df <- cbind(cellranger_se[subjectHits(overlap_se),],macs_se[queryHits(overlap_se),c(1:4)])

write.table(macs_peak_unique,"atac_cuttag_comparison/peak_compare/macs_se_unique.txt",sep = "\t",
            quote = F,row.names = F,col.names = F)
write.table(cellranger_peak_unique,"atac_cuttag_comparison/peak_compare/cellranger_se_unique.txt",sep = "\t",
            quote = F,row.names = F,col.names = F)
write.table(overlap_df,"atac_cuttag_comparison/peak_compare/se_overlap_df.txt",sep = "\t",
            quote = F,row.names = F,col.names = F)


# width histogram
par(mfrow=c(2,2))
cellranger_peak$width <- cellranger_peak$V3-cellranger_peak$V2+1
hist(cellranger_peak$width)
hist(log(cellranger_peak$width))

macs_peak$width <- macs_peak$V3-macs_peak$V2+1
hist(macs_peak$width)
hist(log(macs_peak$width))

cellranger_se$width <- cellranger_se$V3-cellranger_se$V2+1
hist(cellranger_se$width)
hist(log(cellranger_se$width))

macs_se$width <- macs_se$V3-macs_se$V2+1
hist(macs_se$width)
hist(log(macs_se$width))

#------------------------------------------------
# fragments length
#------------------------------------------------
frag_length <- read.table("PBMC/K27ac_PBMC_l1/outs/fragments.tsv",sep="\t",header=F)




