setwd("~/Projects/scCUT_TAG/atac_cuttag_comparison")

atac <- read.table("cluster_compare/pbmc_atac_predicted_id.txt",sep="\t",header=T)
cuttag_cellranger <- read.table("cluster_compare/pbmc_cuttag_predicted_id.txt",sep="\t",header=T)
cuttag_macs <- read.table("cluster_compare/pbmc_cuttag_macs2_predicted_id.txt",sep="\t",header=T)

# CD16+ Monocytes, CD14+ Monocytes, NK cell

atac_bc <- read.table("atac_bc_exmaple",sep = "\t",header=F)
atac_bc_nk <- atac[which(atac$predicted_id == "NK cell"),]
cellranger_bc_nk <- cuttag_cellranger[which(cuttag_cellranger$predicted_id == "NK cell"),]
macs_bc_nk <- cuttag_macs[which(cuttag_macs$predicted_id == "NK cell"),]

cell_group <- unique(c(cuttag_macs$predicted_id,cuttag_cellranger$predicted_id))

final_df <- data.frame()
count_df <- data.frame()
for (c in 1:length(cell_group)){
  print(c)
  cellranger_tmp <- cuttag_cellranger$barcode[which(cuttag_cellranger$predicted_id == cell_group[c])]
  macs_tmp <- cuttag_macs$barcode[which(cuttag_macs$predicted_id == cell_group[c])]
  
  # overlap, and unique list
  common_bc <- intersect(cellranger_tmp,macs_tmp)
  cellranger_unique <- cellranger_tmp[!cellranger_tmp %in%common_bc]
  macs_unique <- macs_tmp[!macs_tmp %in%common_bc]
  
  count_tmp <- data.frame(cell_group=cell_group[c],
                          both=length(common_bc),
                          cellranger=length(cellranger_unique),
                          macs=length(macs_unique))
  tmp_df <- data.frame(barcode=unique(c(cellranger_tmp,macs_tmp)),
                       cell_group=rep(cell_group[c],length(unique(c(cellranger_tmp,macs_tmp)))))
  tmp_df$peak_group <- NA
  for (bc in 1:nrow(tmp_df)){
    if (tmp_df$barcode[bc] %in% common_bc) {
      tmp_df$peak_group[bc] <- "both"
    } else if (tmp_df$barcode[bc] %in% cellranger_unique) {
      tmp_df$peak_group[bc] <- "cellranger"
    } else if (tmp_df$barcode[bc] %in% macs_unique) {
      tmp_df$peak_group[bc] <- "macs2"
    }
  }
  
  final_df <- rbind(final_df,tmp_df)
  count_df <- rbind(count_df,count_tmp)
}

write.table(final_df,"cuttag_bc_cluster_compare.txt",sep="\t",row.names = F,
            quote = F)
write.table(count_df,"cuttag_bc_cluster_compare_count.txt",sep="\t",row.names = F,
            quote = F)
# peaks
cut_cellranger_peak <- read.table("peak_compare/cellranger_peak_unique.txt",sep = "\t",header=F)
overlap_peak <- read.table("peak_compare/overlap_df.txt",sep = "\t",header=F)
macs_peak <- read.table("peak_compare/macs_peak_unique.txt",sep = "\t",header=F)

