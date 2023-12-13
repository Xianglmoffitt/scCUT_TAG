#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# set work directory
setwd("~/Projects/super_enhancer/se_diff_paper/DASE/")
devtools::document()
devtools::load_all()
setwd("~/Projects/super_enhancer/se_data_portal/")
blacklist_df <- read.table("~/Projects/super_enhancer/ENCFF356LFX_blacklist.bed",sep = '\t')

s1 <- args[1]
s2 <- args[2]

# s1 <- cell_combination[1,i]
# s2 <- cell_combination[2,i]

# create dir
out_dir_name <- paste0(s1,"_",s2)
# dir.create(paste0("output/",dir_name))

fc_s1_1_name <- paste0("fc_count/",s1,"_rep1")
fc_s1_2_name <- paste0("fc_count/",s1,"_rep2")
fc_s2_1_name <- paste0("fc_count/",s2,"_rep1")
fc_s2_2_name <- paste0("fc_count/",s2,"_rep2")

# read count table
fc_s1_1 <- read.table(fc_s1_1_name,sep = "\t",header=T)
fc_s1_2 <- read.table(fc_s1_2_name,sep = "\t",header=T)
fc_s2_1 <- read.table(fc_s2_1_name,sep = "\t",header=T)
fc_s2_2 <- read.table(fc_s2_2_name,sep = "\t",header=T)

# change column names
colnames(fc_s1_1)[7] <- "C1_1"
colnames(fc_s1_2)[7] <- "C1_2"
colnames(fc_s2_1)[7] <- "C2_1"
colnames(fc_s2_2)[7] <- "C2_2"

# add enhancer_name
fc_s1_1$enhancer <- paste(fc_s1_1$Chr,fc_s1_1$Start,fc_s1_1$End,sep='_')
fc_s1_2$enhancer <- paste(fc_s1_2$Chr,fc_s1_2$Start,fc_s1_2$End,sep='_')
fc_s2_1$enhancer <- paste(fc_s2_1$Chr,fc_s2_1$Start,fc_s2_1$End,sep='_')
fc_s2_2$enhancer <- paste(fc_s2_2$Chr,fc_s2_2$Start,fc_s2_2$End,sep='_')

# make matrix
fc_s1_new <- merge(fc_s1_1[,c(7,8)],fc_s1_2[,c(7,8)],by="enhancer",all=T)

fc_s2_new <- merge(fc_s2_1[,c(7,8)],fc_s2_2[,c(7,8)],by="enhancer",all=T)

fc_final <- merge(fc_s1_new,fc_s2_new,by="enhancer",all=T)

# replace NA with 0
fc_final[is.na(fc_final)] <- 0

# save count table
write.table(fc_final,file=paste0("output/",out_dir_name,"/","count_table.txt"),
            sep="\t",row.names = F,quote = F)

#-----------------------------------------------------------------------------
# run DASE with count matrix
#-----------------------------------------------------------------------------
# enhancer
en_s1_1 <- read.table(paste0("enhancer_bed/",s1,"_rep1_peaks.narrowPeak"),sep='\t',header =F)
en_s1_2 <- read.table(paste0("enhancer_bed/",s1,"_rep2_peaks.narrowPeak"),sep='\t',header =F)
en_s2_1 <- read.table(paste0("enhancer_bed/",s2,"_rep1_peaks.narrowPeak"),sep='\t',header =F)
en_s2_2 <- read.table(paste0("enhancer_bed/",s2,"_rep2_peaks.narrowPeak"),sep='\t',header =F)
pool_enhancer_df <- unique(rbind(en_s1_1,en_s1_2,en_s2_1,en_s2_2))

# SE
se_s1_1 <- read.table(paste0("se_bed/",s1,"_rep1_peaks_Gateway_SuperEnhancers.bed"),sep='\t',header =F)
se_s1_2 <- read.table(paste0("se_bed/",s1,"_rep2_peaks_Gateway_SuperEnhancers.bed"),sep='\t',header =F)
se_s2_1 <- read.table(paste0("se_bed/",s2,"_rep1_peaks_Gateway_SuperEnhancers.bed"),sep='\t',header =F)
se_s2_2 <- read.table(paste0("se_bed/",s2,"_rep2_peaks_Gateway_SuperEnhancers.bed"),sep='\t',header =F)
pool_se_df <- rbind(se_s1_1,se_s1_2,se_s2_1,se_s2_2)

# run DASE
dase_out <- DASE(se_in = pool_se_df,
                 e_in = pool_enhancer_df,
                 enhancer_count_table = fc_final,
                 c1_n=2,
                 c2_n=2,
                 bl_file = blacklist_df)
#-----------------------------------------------------------------------------
# save results
#-----------------------------------------------------------------------------
saveRDS(dase_out, file = paste0("output/",out_dir_name,"/",out_dir_name,"_dase_out.rds"))

write.table(dase_out$se_category,file = paste0("output/",out_dir_name,"/se_category.txt"),sep='\t',
            row.names = FALSE,quote = FALSE)
write.table(dase_out$ce_fit,file = paste0("output/",out_dir_name,"/ce_fit.txt"),sep='\t',
            row.names = FALSE,quote = FALSE)
write.table(dase_out$cutoff,file = paste0("output/",out_dir_name,"/cutoff.txt"),sep='\t',
            row.names = FALSE,quote = FALSE)

pdf(file = paste0("output/",out_dir_name,"/boxplot.pdf"))
dase_out$boxplot
dev.off()
pdf(file = paste0("output/",out_dir_name,"/density.pdf"))
dase_out$density_plot
dev.off()

pdf(paste0("output/",out_dir_name,"/MAplot.pdf"))
plotMA(dase_out$lfc_shrink)
dev.off()


