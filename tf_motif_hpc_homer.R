#!/usr/bin/env Rscript

library(GenomicRanges)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)

# rose peak  table
#se_profile_out<- readRDS(args[1])
peak_bed <- read.table(args[1],header=F,sep="\t")
peak_bed$peak <- paste(peak_bed$V2,peak_bed$V3,peak_bed$V4,sep="_")
target_chr <- c("chr1","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chr10","chr11","chr12",
                "chr13","chr14","chr15","chr16","chr17","chr18",
                "chr19","chr20","chr21","chr22","chrX")

# only get target chr
peak_bed <- peak_bed[which(peak_bed$V2 %in% target_chr),]

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
  seqnames=peak_bed$V2,
  ranges=IRanges(peak_bed$V3,peak_bed$V4,names=peak_bed$peak)
)

# extract sequence
hg38 <- BSgenome.Hsapiens.UCSC.hg38
seq <- getSeq(hg38,gr_e)

# get tfbs for each enhancer
# parallel
registerDoParallel(5)

#start <- proc.time()
se_tf <- foreach (i= 1:length(seq), .combine = rbind) %dopar% {
  search_temp<-searchSeq(pwm_jaspar,
                         seq[i],
                         min.score = 0.9,
                         strand="*")
  ma_count_temp<- sapply(search_temp, function(x) length(x))
  ma_count <- which(ma_count_temp !=0)
  se_tf_temp <- data.frame(e_merged_name = rep(names(seq)[i],length(ma_count)),
                           ma_id = gsub(paste0("\\.",names(seq)[i]),"",names(ma_count)))
  
  se_tf_temp
}

write.table(se_tf,file = args[2],sep='\t',
            row.names = FALSE,quote = FALSE)

#proc.time()-start
