#!/usr/bin/env Rscript

setwd("~/Projects/scCUT_TAG/cut_chip_comp/")
library(GenomicRanges)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)

# rose peak  table
peak_bed <- read.table("macs_out/H1_chip_1_peaks.narrowPeak",header=F,sep="\t")
#peak_bed <- read.table(args[1],header=F,sep="\t")


target_chr <- c("chr1","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18",
                "chr19","chr20","chr21","chr22","chrX")

# only get target chr
peak_bed <- peak_bed[which(peak_bed$V1 %in% target_chr),]

peak_bed$peak <- paste(peak_bed$V1,peak_bed$V2,peak_bed$V3,sep="_")
# create jaspar TF database
library(JASPAR2018)
pfm_jaspar <- getMatrixSet(JASPAR2018,
                           opts=list(
                             collection = 'CORE',
                             species='9606'
                           ))
pwm_jaspar <- toPWM(pfm_jaspar,type="log2probratio")

# create enhancer range
# enhancer_fit$g_start <- enhancer_fit$start+enhancer_fit$width/2-500
# enhancer_fit$g_end <- enhancer_fit$start+enhancer_fit$width/2+500

# extract significant enhancer first to reduce computational time
peak_bed <- head(peak_bed)
gr_e <- GRanges(
  seqnames=peak_bed$V1,
  ranges=IRanges(peak_bed$V2,peak_bed$V3,names=peak_bed$peak)
)

# extract sequence
hg38 <- BSgenome.Hsapiens.UCSC.hg38
seqlist <- getSeq(hg38,gr_e)

search_temp<-searchSeq(pwm_jaspar,
                       seqlist,
                      min.score = 0.9,
                      strand="*")
se_gff <- writeGFF3(search_temp)

write.table(se_gff,file = args[2],sep='\t',
            row.names = TRUE,quote = FALSE)

# se_tf <- data.frame()
# for (i in seq(1,length(seqlist),200)) {
#   print(i)
#   if (i ==10601) {
#     search_temp<-searchSeq(pwm_jaspar,
#                            seqlist[i:length(seqlist)],
#                            min.score = 0.9,
#                            strand="*")
#   } else {
#     search_temp<-searchSeq(pwm_jaspar,
#                            seqlist[i:(i+200)],
#                            min.score = 0.9,
#                            strand="*")
#   }
# 
#   set_tmp <- writeGFF3(search_temp)
#   set_tmp$ma_id <- gsub(paste0("\\.",names(seq)[i]),"",row.names(set_tmp))
#   se_tf <- rbind(se_tf,set_tmp)
# }

# write.table(se_tf,file = args[2],sep='\t',
#             row.names = FALSE,quote = FALSE)


# get tfbs for each enhancer
# parallel
# registerDoParallel(5)
# 
# #start <- proc.time()
# se_tf <- foreach (i= 1:length(seq), .combine = rbind) %dopar% {
#     search_temp<-searchSeq(pwm_jaspar,
#                          seq[i],
#                          min.score = 0.9,
#                          strand="*")
#     ma_count_temp<- sapply(search_temp, function(x) length(x))
#     ma_count <- which(ma_count_temp !=0)
#     se_tf_temp <- data.frame(e_merged_name = rep(names(seq)[i],length(ma_count)),
#                            ma_id = gsub(paste0("\\.",names(seq)[i]),"",names(ma_count)))
#     se_tf_temp
# }



#proc.time()-start
