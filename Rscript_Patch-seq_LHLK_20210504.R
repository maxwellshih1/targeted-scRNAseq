if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("RnaSeqSampleSize", force = TRUE)
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("apeglm", force = TRUE)
BiocManager::install("glmGamPoi")
install.packages("tidyverse")
install.packages("data.table")
install.packages("RColorBrewer")
install.packages("viridis")
# install.packages("ggplot2")
# install.packages("tidyr")
install.packages("pheatmap")
BiocManager::install("edgeR")
# BiocManager::install("statmod")
# library(RnaSeqSampleSize)
# library(RnaSeqSampleSizeData)
library("BiocParallel")
library(DESeq2)
library(apeglm)
library(glmGamPoi)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(data.table)
# library(edgeR)
# library(statmod)
# Remove NA by zero using Matt Dowle's function 
#ref: https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table/7249454#7249454
f_dowle2 = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=0]
}
# browseVignettes("RnaSeqSampleSize")



getwd()
setwd("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev")
# sc-qRT-PCR analysis ####
fread("./LHLK_sc-qRT-PCR_data_20171215/Mx_analysis-Table 1.csv", col.names = c("sample_id","Fly_id","Hemi","gender","Act5c_r1","Act5c_r2", "Act5c mean","Lk_r1", "Lk_r2","Lk mean", "coGFP_r1","coGFP_r2", "coGFP mean","deltaCт Mean (Lk-Act5c)", "deltaCт Mean (coGFP-Act5c)")) -> qpcr_dt
names(qpcr_dt) 
qpcr_dt[,`:=`(sample_id=as.factor(sample_id),
              Fly_id=as.factor(Fly_id),
              Hemi=as.factor(Hemi),
              gender=as.factor(gender),
              Lk_r1=as.numeric(Lk_r1),
              Lk_r2=as.numeric(Lk_r2),
              coGFP_r1=as.numeric(coGFP_r1))][]
qpcr_dt[,.SD,.SDcol=c(1:3,5:6,8:9,11:12)] -> ht_map_qpcr_dt
t(data.frame(ht_map_qpcr_dt[,.SD,.SDcol=c(grep("r1$|r2$",names(ht_map_qpcr_dt),perl = T))],row.names = ht_map_qpcr_dt$sample_id)) -> ht_map_df
data.frame(ht_map_qpcr_dt[,.SD,.SDcol=c(2:3)],row.names = ht_map_qpcr_dt$sample_id) -> ht_map_colData
as.matrix(ht_map_df)
ht_map_clr = list(Fly_id=c(`1`="#7fc97f",`3`="#beaed4",`4`="#fdc086",`5`="#ffff99",`6`="#386cb0",`7`="#f0027f",`8`="#bf5b17"), Hemi=c(L="#DFDF00", R="#009F7F" ) )
# ph.clr = list( Cell_type=c(MBON.g3bp1="#FF00FF",LHLK="#1FDFFF"), State=c(Fed="#DFDF00", Starved="#009F7F" ) )

curBreaks <- seq(min(ht_map_df[!is.na(ht_map_df)]), max(ht_map_df[!is.na(ht_map_df)]), length.out=101)

# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/single-cell_qRT-PCR/col_KeeneLab/ht-map_10lhlk_qRT-PCR_01JUN2024.pdf",height=4)
pheatmap(ht_map_df[,c(1,3,5,7,9,2,4,6,8,10)], color=viridis_pal(begin=0.3,direction =-1, option ="A")(101),cluster_rows=F, show_rownames=T, cluster_cols=F,  legend = T, annotation_col = ht_map_colData, annotation_colors=ht_map_clr, breaks = curBreaks,
         main = "", display_numbers = F, number_color = "grey90",angle_col = 0,
         fontsize = 12, cellwidth = NA, treeheight_row = NA, treeheight_col = 50, na_col="#000000")
dev.off()





# Import the joined count tab ####
## Lk.1
lk1.counts <- fread("v4_lk1.counts") 

grep("^v4_lk1_",list.files(), value = T, perl = T) -> lk1_id
gsub("^.*([0-9]{3}).*$","lk1_\\1",lk1_id) -> lk1_id
names(lk1.counts) <- c('gene_id',lk1_id)
# is.na(lk1.counts) %>%  sum()
lk1.counts[!grep('__',gene_id),] -> lk1.counts

lk1.counts[gene_id=="Lk",]
colSums(lk1.counts[,2:20]) %>% as.numeric()
# fwrite(lk1.counts, "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/lk1_counts.tab")

## Lk.3
lk3.counts <- fread("v4_lk3.counts")

grep("^v4_lk3_",list.files(), value = T, perl = T) -> lk3_id
gsub("v4_","",lk3_id) -> lk3_id
gsub("_R2","",lk3_id) -> lk3_id
gsub(".dir.counts","",lk3_id) -> lk3_id
names(lk3.counts) <- c('gene_id',lk3_id)
# is.na(lk3.counts) %>%  sum()
lk3.counts[grep('__',gene_id),]
lk3.counts[!grep('__',gene_id),] -> lk3.counts
lk3.counts[gene_id=="Lk",]
colSums(lk3.counts[,2:19]) %>% as.numeric()
# fwrite(lk3.counts, "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/lk3_counts.tab")

## Lk.4
lk4.counts <- fread("v4_lk4.counts")

grep("^v4_lk4_",list.files(), value = T, perl = T) -> lk4_id
gsub("v4_","",lk4_id) -> lk4_id
gsub("_R2","",lk4_id) -> lk4_id
gsub(".dir.counts","",lk4_id) -> lk4_id
names(lk4.counts) <- c('gene_id',lk4_id)
# is.na(lk4.counts) %>%  sum()

lk4.counts[!grep('__',gene_id),] -> lk4.counts
lk4.counts[gene_id=="Lk",]
colSums(lk4.counts[,2:29]) %>% as.numeric()
# fwrite(lk4.counts, "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/lk4_counts.tab")

# Normalization ####
CountToCpm <- function(Count)
{
  exp(log(Count) - log(sum(Count)) + log(1e6))
}
data.table(gene_id = lk1.counts$gene_id,
           lk1.counts[,lapply(.SD,CountToCpm),.SDcols=2:ncol(lk1.counts)]) -> lk1_tpm
colSums(lk1_tpm[,-1])
# fwrite(lk1_tpm, "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/lk1_tpm.tab")
data.table(gene_id = lk3.counts$gene_id,
           lk3.counts[,lapply(.SD,CountToCpm),.SDcols=2:ncol(lk3.counts)]) -> lk3_tpm
colSums(lk3_tpm[,-1])
# fwrite(lk3_tpm, "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/lk3_tpm.tab")
data.table(gene_id = lk4.counts$gene_id,
           lk4.counts[,lapply(.SD,CountToCpm),.SDcols=2:ncol(lk4.counts)]) -> lk4_tpm
colSums(lk4_tpm[,-1])
# fwrite(lk4_tpm, "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/lk4_tpm.tab")

#### Calculate the ERCC ####
e.t.ratio <- c( colSums(lk1.counts[grep("ERCC",gene_id),2:20])/colSums(lk1.counts[!grep("ERCC",gene_id),2:20]), 
                colSums(lk3.counts[grep("ERCC",gene_id),2:19])/colSums(lk3.counts[!grep("ERCC",gene_id),2:19]),
                colSums(lk4.counts[grep("ERCC",gene_id),2:29])/colSums(lk4.counts[!grep("ERCC",gene_id),2:29]))
str(e.t.ratio)
as.numeric(e.t.ratio)
as.data.frame(e.t.ratio) -> df
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_bar_all.lk_SpikeInratio.20210504.pdf",
#     width = 17)
# ggplot(data = df, aes(x = rownames(df), y=e.t.ratio)) +
#   geom_bar(stat = "identity", width = 0.6) + #scale_x_discrete(limits=rownames(colData)) +
#   theme_gray(base_size = 16) + theme(axis.text.x = element_text(angle=45, hjust = 1,
#                                                                 vjust = 1)) +
#   labs(title="(0.09 cutoff)", x ="Cell id", y = "Spike-in read ratio") +
#   geom_hline(yintercept=0.09, linetype="dashed", color = "red")
# dev.off()

## Spike-in linearity ####
fread("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/seq.data/ERCC92_conc.txt", header = T, sep = "\t") -> ERCC.conc
names(ERCC.conc) <- c("ERCC.ID", "attoMoles.o.micron")
head(ERCC.conc)
tail(ERCC.conc)
mutate(ERCC.conc,
       molecule.o.micron = attoMoles.o.micron * 6.022140 * 10^23/10^18,
       spike.in.lk = molecule.o.micron*0.4/(4*10^5) ) -> ERCC.conc
str(ERCC.conc)
head(ERCC.conc)
# ERCC.conc[ERCC.ID=="ERCC-00017",]
spike.in.corr <- NA
lk4_tpm[ERCC.conc, on=c("gene_id==ERCC.ID")] -> spike.in.corr
lk3_tpm[spike.in.corr, on=c("gene_id==gene_id")] -> spike.in.corr
lk1_tpm[spike.in.corr, on=c("gene_id==gene_id")] -> spike.in.corr
tail(spike.in.corr)
duplicated(spike.in.corr$gene_id)
length(spike.in.corr$gene_id)
# Replacde the NA with 0
f_dowle2(spike.in.corr)
anyNA(spike.in.corr)
# In CEL-Seq paper (Hashimshony et al., 2012), log scale was used to calculate the correlation
spike.in.corr
f_mx1 <- function(tpm)
{
  log10(tpm+1)
}
data.table(ERCC_ID=spike.in.corr$gene_id,
           spike.in.corr[,lapply(.SD,f_mx1),.SDcols=2:ncol(spike.in.corr)]) -> l.spike.in.corr
l.spike.in.corr[,]
nrow(l.spike.in.corr)
spikein.lin <-c()
data.frame(l.spike.in.corr[,2:ncol(l.spike.in.corr)], row.names = l.spike.in.corr$ERCC_ID) -> l.spike.cor.df
for (i in 1:length(c(lk1_id, lk3_id, lk4_id)) ) {
  outfile <- paste("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_scatter_ERCC.corr_", names(l.spike.cor.df)[i], ".pdf", sep = "" )
  gg <- ggplot(l.spike.cor.df, aes(x = spike.in.lk, y = l.spike.cor.df[,i] ) )
  gg.Ecc=round(cor(l.spike.cor.df[,"spike.in.lk"] , l.spike.cor.df[,i]),digits = 3)
  gg.Ecc-> spikein.lin[i]
  # pdf(outfile)
  # print(gg + geom_point(size = 2) + geom_smooth(method = lm, se = F) + theme_gray(base_size = 20) + 
  #         labs(title=paste("Linearity of spike-in for ",names(l.spike.cor.df)[i], ", R=", gg.Ecc, sep = ""), 
  #              x ="Spike-in introduced", y = "Spike-in detected") )
  # dev.off()
}
spikein.lin

# Remove the ERCC 
lk1.counts
lk1.counts[!grepl('ERCC', lk1.counts$gene_id ),] -> lk1.counts
lk1_tpm[!grepl('ERCC', lk1_tpm$gene_id),] -> lk1_tpm
lk3.counts[!grepl('ERCC', lk3.counts$gene_id ),] -> lk3.counts
lk3_tpm[!grepl('ERCC', lk3_tpm$gene_id),] -> lk3_tpm
lk4.counts[!grepl('ERCC', lk4.counts$gene_id ),] -> lk4.counts
lk4_tpm[!grepl('ERCC', lk4_tpm$gene_id),] -> lk4_tpm
# Count the detected genes
ncol(lk1.counts)
ncol(lk3.counts)
ncol(lk4.counts)
d.g.num <- c(colSums(lk1.counts[,2:20]>0), colSums(lk3.counts[,2:19]>0), colSums(lk4.counts[,2:29]>0) )
as.vector(d.g.num)

## Read-in the bootstramped colData tab ####
fread("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/colData.lk.tsv") -> colData.dt

data.frame(colData.dt, row.names = c(lk1_id, lk3_id, lk4_id)) -> colData.df
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000 & 
             colData.df$Cell_id!=45, ] -> g.colData.df # In DPM project, 0.035 and 0.8 as thresholds
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000 & 
             colData.df$Cell_id!=45, "Train_protocol"] %>% table()
## Heat mapping ####
lk1_tpm
lk1_tpm[lk3_tpm, on=c("gene_id==gene_id")] -> lk_tpm
lk_tpm[lk4_tpm, on=c("gene_id==gene_id")] -> lk_tpm
f_dowle2(lk_tpm)
data.frame(lk_tpm[,2:(length(c(lk1_id, lk3_id, lk4_id))+1)],row.names = lk_tpm$gene_id) -> lk_tpm.df

colData.df[,c("Cell_type", "Train_protocol")] -> ph.colData
ph.colData["lk1_041",1]<- c("MBON.g3bp1")
names(ph.colData) <- c('Cell_type',"State")
ph.clr = list( Cell_type=c(MBON.g3bp1="#FF00FF",LHLK="#1FDFFF"), State=c(Fed="#DFDF00", Starved="#009F7F" ) )

# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.r-cluster.lk.1-4.20210504.pdf")
# pheatmap(log10(lk_tpm.df+1), cluster_rows=T, show_rownames=F, cluster_cols=F,  legend = T, annotation_col = ph.colData, annotation_colors = ph.clr,
#          main = "Transcriptome profiles of Lk.1 & Lk.3 & Lk.4 (log TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
#          fontsize = 7.5, cellwidth = NA, treeheight_row = 0, treeheight_col = 50)
# dev.off()
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.rc-cluster.lk.1-4.20210504.pdf")
# pheatmap(log10(lk_tpm.df+1), cluster_rows=T, show_rownames=F, cluster_cols=T,  legend = T, annotation_col = ph.colData, annotation_colors = ph.clr,
#          main = "Transcriptome profiles of Lk.1 & Lk.3 & Lk.4 (log TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
#          fontsize = 7.5, cellwidth = NA, treeheight_row = 0, treeheight_col = 50)
# dev.off()
# 
# data.table(colData.df)[,unique(Gender)]
# filter(colData.df, Cell_type == 'LHLK') %>% select(Gender, Train_protocol) -> ph.colData
# names(ph.colData) <- c('Gender',"State")
# ph.clr = list( Gender=c(Female="#FF00FF",Male="#1FDFFF"), State=c(Fed="#DFDF00", Starved="#009F7F" ) )
# 
# head(lk_tpm.df)
# select(lk_tpm.df, rownames(ph.colData)) -> LHLK_tpm.df
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap_Gender-State_lk.1-4.20240317.pdf")
# pheatmap(log10(LHLK_tpm.df+1), cluster_rows=T, show_rownames=F, cluster_cols=F,  legend = T, annotation_col = ph.colData, annotation_colors = ph.clr,
#          main = "Transcriptome profiles of LHLK (log TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
#          fontsize = 7.5, cellwidth = NA, treeheight_row = 0, treeheight_col = 50)
# dev.off()




# Quality control ####

# Lk.1 + Lk.3 + Lk.4 with filter of 
#Cell_type=="LHLK" & ERCCoverTranscript_ratio<0.09 & Spike.in.linearity>0.80 & Assigned.reads >250000
## lk.1_045 has no signal for coGFP!! # lk.4_029 is from a male fly

colData.df[colData.df$Cell_type=="LHLK"&colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000, ] %>% nrow()
colData.df[colData.df$Cell_type=="LHLK"&colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000, ] %>% rownames() -> g_42LHLK_id
# 42 LHLK \ 24 Fed & 18 Starved
# lk_tpm.df %>% select(all_of(g_42LHLK_id)) -> lk_tpm.df
# lk_tpm.dt %>% select(all_of(g_42LHLK_id)) -> lk_tpm.dt
# lk_tpm %>% select(all_of(g_42LHLK_id)) -> lk_tpm
# lk_counts %>% select(all_of(g_42LHLK_id)) -> lk_counts

# ## Sanity check heatmap per Josh's request ####
# # bigger font
# # "coGFP","Lk","trsn","elav","repo","roX2" # remove trsn, roX2
# data.frame(id = c(paste("Lk.",1:24, sep = "")),row.names =  g.lk.fed.id) -> ph.col.id
# sanity.c <- c("coGFP","Lk","elav","repo")
# g.lk.fed_tpm.df[c(sanity.c),] -> sanity.g.lk
# # round(log10(sanity.g.lk+1), digits = 1) -> m
# curBreaks <- seq(0, 3, length.out=101)
# names(sanity.g.lk) == rownames(ph.col.id)
# ncol(sanity.g.lk)
# # pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.lk_24.fed.gLHLKs_sanity.log.20210504.pdf", width = 11, height = 2.7)
# # pheatmap(log10(sanity.g.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col= NA, labels_col = c(1:24), breaks = curBreaks,
# #          legend = T, angle_col = 90,
# #          main = "Marker Gene Expression", display_numbers = F, #number_color = "grey10", number_format = "%.1f",
# #          #          fontsize_col = 20,
# #          fontsize = 17, treeheight_row = 0, cellwidth = NA)
# # dev.off()

# #ERCCoverTranscript_ratio
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_bar_all.lk_SpikeInR.20210504.pdf", width = 22)
# ggplot(data = colData.df, aes(x = rownames(colData.df), y=ERCCoverTranscript_ratio)) + geom_bar(stat = "identity", width = 0.6) +
#   scale_x_discrete(limits=rownames(colData.df)) + theme_gray(base_size = 18) +
#   theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) +
#   labs(title="The higher, the more RNA degraded (cutoff = 0.09)", x ="Cell_id", y = "Spike-in to transcript ratio") +
#   geom_hline(yintercept=0.09, linetype="dashed", color = "red")
# dev.off()

# #colData.df$Spike.in.linearity
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_bar_all.lk_SpikeInLinearity.20210504.pdf", width = 22)
# ggplot(data = colData.df, aes(x = rownames(colData.df), y=Spike.in.linearity)) + geom_bar(stat = "identity", width = 0.6) +
#   scale_x_discrete(limits=rownames(colData.df)) + theme_gray(base_size = 18) +
#   theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) +
#   labs(title="Correlation between numbers of spike-in molecules and quantification (cutoff = 0.8)", x ="Cell_id", y = "Spike-in linearity") +
#   geom_hline(yintercept=0.8, linetype="dashed", color = "red")
# dev.off()
colData.df$Assigned.reads
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_bar.log_all.lk_AssignedRead.20210504.pdf", width = 22)
# ggplot(data = colData.df, aes(x = rownames(colData.df), y=Assigned.reads)) + geom_bar(stat = "identity", width = 0.6) +
#   scale_x_discrete(limits=rownames(colData.df)) + scale_y_log10(breaks = c(10^(0:7))) + coord_cartesian(xlim = NULL,ylim = c(10^4, 10^6.5)) +
#   theme_gray(base_size = 18) + theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) +
#   labs(title="250k is sufficient for good quantification in CEL-Seq1 paper", x ="Cell_id", y = "Assigned read number") +
#   geom_hline(yintercept=250000, linetype="dashed", color = "red")
# dev.off()
#colData.df$detected.gene.count
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_bar_all.lk_detected.gene.count.20210504.pdf", width = 22)
# ggplot(data = colData.df, aes(x = rownames(colData.df), y=detected.gene.count)) + geom_bar(stat = "identity", width = 0.6) +
#   scale_x_discrete(limits=rownames(colData.df)) + #scale_y_continuous(labels = c('0','20k','40k','60k','80k')) +
#   theme_gray(base_size = 18) + theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) +
#   labs(title="", x ="Cell_id", y = "Detected gene count")
# dev.off()
#colData.df$algn.rate
parse_number(colData.df$algn.rate)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_bar_all.lk_algn.rate.20210504.pdf", width = 22)
# ggplot(data = colData.df, aes(x = rownames(colData.df), y=parse_number(colData.df$algn.rate))) + geom_bar(stat = "identity", width = 0.6) +
#   scale_x_discrete(limits=rownames(colData.df)) + #scale_y_continuous(labels = c('0','20k','40k','60k','80k')) +
#   theme_gray(base_size = 18) + theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) +
#   labs(title="", x ="Cell_id", y = "Uniquely aligned rate (%)")
# dev.off()
#colData.df$Counted.reads
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_bar_all.lk_CountedRead.20210504.pdf", width = 22)
# ggplot(data = colData.df, aes(x = rownames(colData.df), y=uniq.reads)) + geom_bar(stat = "identity", width = 0.6) +
#   scale_x_discrete(limits=rownames(colData.df)) + scale_y_continuous(labels = c('0','50k','100k','150k')) +
#   theme_gray(base_size = 18) + theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) +
#   labs(title="median 162k for SCRSII-DPM dataset (cutoff = 20,000)", x ="Cell_id", y = "Unique read count") +
#   geom_hline(yintercept=20000, linetype="dashed", color = "red")
# dev.off()

behavior_gene <- c("inc", "regucalin", "deltaCOP", "Sec61gamma", "Sec61beta", "inc", "fit", "Sodh1" )
# 44 lhlk? ####
nrow(colData.df)

# # 24 Fed vs 18 Starved ####
# nrow(colData.df[g_42LHLK_id,])
# colData.df[g_42LHLK_id,]-> g_42LHLK_colData.df
# 
# # State=c(Fed="#DFDF00", Starved="#009F7F" )
# g_42LHLK_colData.df %>% select(Cell_id ,Raw_read_num, uniq.reads, Train_protocol, detected.gene.count) -> est.depth
# class(est.depth$detected.gene.count)
# as.numeric(est.depth$uniq.reads) -> est.depth$uniq.reads
# # pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_scatter_42LHLKs_detected.gene-log.uni.r.20210504.pdf",
# #     width = 15, height = 8)
# # sc <- ggplot(est.depth, aes(x = uniq.reads, y = detected.gene.count) )
# # sc + geom_point(aes(color = Train_protocol),size= 2.5) +
# #   scale_color_manual(breaks = c("Fed", "Starved"),values=c("#DFDF00","#009F7F")) +
# #   geom_smooth(method = 'lm', se = F, formula = y ~ log(x), color = 'light blue',alpha = 0.2 ) +
# #   geom_text(aes(label = Cell_id ), size = 6, vjust = 0, hjust = 0, alpha = 0.65,check_overlap = T ) +
# #   theme_gray(base_size = 20) + scale_x_log10() +
# #   labs(title=paste("24 Fed LHLKs & 18 Starved LHLKs\nmedian detected gene count ",
# #                    median(as.numeric(est.depth$detected.gene.count)),"\nmedian unique read count ",
# #                    median(as.numeric(est.depth$uniq.reads)), sep = ''),
# #        x ="log10(unique read count)", y = "Detected gene count" ) +
# #   geom_hline(yintercept= median(as.numeric(est.depth$detected.gene.count)), linetype="dashed", color = "black", alpha = 0.5) +
# #   geom_vline(xintercept= median(as.numeric(est.depth$uniq.reads)), linetype="dashed", color = "black", alpha = 0.5)
# # dev.off()
# ###
# 
# ## DESeq2 for 18 starved vs 24 fed ####
# ncol(count.df)
# count.df[,g_42LHLK_id] -> g_42LHLK_count.df
# ncol(g_42LHLK_count.df)
# g_42LHLK_colData.df
# str(g_42LHLK_colData.df)
# class(g_42LHLK_colData.df$Library)
# as.factor(g_42LHLK_colData.df[,"Train_protocol"]) -> g_42LHLK_colData.df$Train_protocol
# 
# ### (obsolete) Train_protocl only DESeq ####
# # decided to go Library + Train_protocol
# dds <- DESeqDataSetFromMatrix(countData = g_42LHLK_count.df,
#                               colData = g_42LHLK_colData.df,
#                               design = ~ Train_protocol)
# # First, we filter out the lowly expressed genes, by removing those genes that 
# #   do not have at least 1 reads in at least half of the samples.
# dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
# # filter <- rowSums(assay(dds)>1)>5
# # table(filter) # FALSE 10281, TRUE 6830
# # Identify the 2000 most variable genes for DE analysis will miss some potential DE genes!! 
# assay(dds.preF) %>% log1p %>% rowVars(useNames = T) -> vars
# # str(vars)
# # names(vars) <- rownames(dds.preF)
# vars <- sort(vars, decreasing = TRUE) 
# head(vars)
# dds.preF <- dds.preF[names(vars)[1:2000],]
# 
# dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref='Fed')
# dds.preF <- DESeq(dds.preF)
# resultsNames(dds.preF)
# res <- results(dds.preF, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
# res
# summary(res)
# res.Ordered <- res[order(res$padj),]
# res.Ordered$padj <- ifelse(is.na(res.Ordered$padj), 1, res.Ordered$padj)
# # write.csv(res.Ordered[Lib_res.Ordered$padj<=0.050,],"/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/DE_42lhlk_state_18MAR2024.csv")
# g42lhlk.gene = rownames(res.Ordered[res.Ordered$padj<=0.050,])
# g42lhlk.gene
# 
# plotMA(res, ylim=c(-2,2))
# 
# 
# ### Library + Train_protocl only DESeq ####
# dds <- DESeqDataSetFromMatrix(countData = g_42LHLK_count.df, 
#                               colData = g_42LHLK_colData.df, 
#                               design = ~ Library + Train_protocol)
# # First, we filter out the lowly expressed genes, by removing those genes that 
# #   do not have at least 1 reads in at least half of the samples.
# dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
# # Identify the 2000 most variable genes for DE analysis will miss some potential DE genes!! 
# assay(dds.preF) %>% log1p %>% rowVars(useNames = T) -> vars
# # str(vars)
# # names(vars) <- rownames(dds.preF)
# vars <- sort(vars, decreasing = TRUE) 
# head(vars)
# dds.preF <- dds.preF[names(vars)[1:2000],]
# # keep <- rowSums(counts(dds)) >= 20 # stick to the original prefiltering mean 2024.03.23
# # dds <- dds[keep,]
# dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref='Fed')
# dds.preF <- DESeq(dds.preF, test="LRT", reduced = ~Library, useT=T, minmu=1e-6,
#              minReplicatesForReplace = Inf) # use single-cell recommendation settings
# resultsNames(dds.preF)
# Lib_res <- results(dds.preF, contrast = c("Train_protocol","Starved","Fed")) #alpha=0.05
# Lib_res
# summary(Lib_res)
# Lib_res.Ordered <- Lib_res[order(Lib_res$padj),]
# Lib_res.Ordered$log2FoldChange %>% summary()
# # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# # -2.292718 -0.238993 -0.006987 -0.012374  0.213073  2.644323 
# Lib_res.Ordered$padj <- ifelse(is.na(Lib_res.Ordered$padj), 1, Lib_res.Ordered$padj) # to remove the error message
# Lib_res.Ordered['inc',]
# # Error: logical subscript contains NAs
# # write.csv(Lib_res.Ordered, "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/sc.2kVars.DE_42lhlk_lib-state_23MAR2024.csv")
# Lib_g42lhlk.gene = rownames(Lib_res.Ordered[Lib_res.Ordered$padj<=0.10,])
# Lib_g42lhlk.gene
# 
# DESeq2::plotMA(Lib_res, ylim=c(-3,3), ylab="log2 fold change", main = "alpha = 0.1")
# dev.off()
# #### heat map ####
# colData.df[c(g_42LHLK_id),c("Gender","Library", "Train_protocol")] -> ph.colData
# names(ph.colData)[3] <- c("State")
# rownames(ph.colData[ph.colData$State=='Fed',]) -> fed_24lhlk_id
# setdiff(g_42LHLK_id, fed_24lhlk_id) -> starved_18lhlk_id
# # brute force to find the col_order
# lk_tpm.df[c(Lib_g42lhlk.gene),c(fed_24lhlk_id)] -> de_24.lk
# log10(de_24.lk+1)
# curBreaks <- seq(0, 3, length.out=101)
# 
# ph.clr = list( Gender=c(Female="#FF00FF",Male="#1FDFFF"),
#                Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
#                State=c(Fed="#DFDF00", Starved="#009F7F" ) )
# 
# pheatmap(log10(de_24.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
#          color=viridis_pal(option = "A")(101) ,fontsize = 7.5, treeheight_row = 0, cellwidth = NA)
# dev.off()
# 
# lk_tpm.df[c(Lib_g42lhlk.gene),c(starved_18lhlk_id)] -> de_18.lk
# log10(de_18.lk+1)
# curBreaks <- seq(0, 3, length.out=101)
# 
# ph.clr = list( Gender=c(Female="#FF00FF",Male="#1FDFFF"),
#                Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
#                State=c(Fed="#DFDF00", Starved="#009F7F" ) )
# 
# pheatmap(log10(de_18.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
#          color=viridis_pal(option = "A")(101) ,fontsize = 7.5, treeheight_row = 0, cellwidth = NA)
# dev.off()
# cluster_order_42lhlk_id = gsub("^","lk",c('3_096','1_058','3_089','1_055','1_057',
#                                           '1_054','3_088','4_099','4_111','3_078',
#                                           '4_115','4_110','4_118','4_098','4_116',
#                                           '4_113','4_117','4_119','4_122','4_112',
#                                           '4_121','4_029','1_046','4_114',
#                                           '1_043','1_052','1_049','3_100','1_053',
#                                           '1_044','1_050','3_086','4_105','4_028',
#                                           '4_123','4_103','4_107','4_104','3_101',
#                                           '4_102','3_092','4_108'))
# lk_tpm.df[c(Lib_g42lhlk.gene),c(cluster_order_42lhlk_id)] -> de_42.lk
# log10(de_42.lk+1)
# max(log10(de_42.lk+1)) #3.3527
# min(log10(de_42.lk+1)) #0
# curBreaks <- seq(0, 3, length.out=101)
# names(ph.colData)[3] <- c("State")
# ph.clr = list( Gender=c(Female="#FF00FF",Male="#1FDFFF"),
#                Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
#                State=c(Fed="#DFDF00", Starved="#009F7F" ) )
# 
# # pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap_42lhlk_sc.2kVars.DESeq-24genes_23MAR2024.pdf", width = 17)
# # pheatmap(log10(de_42.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T, legend_breaks = c(0:3), 
# #          main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
# #          color=viridis_pal(option = "A")(101) ,fontsize = 15, treeheight_row = 0, cellwidth = NA)
# # dev.off()

# 23 feb vs 18 starved ####
### Time_anesthetized.h ####
Time_anesthetized.h = c(1.00,1.00,2.75,3.60,3.60,3.60,0.75,2.00,2.00,3.85,3.85,3.60,2.70,4.10,4.40,6.60,2.25,3.75,4.40,2.50,2.65,2.65,4.40,1.40,1.70,1.70,2.80,2.80,4.30,4.55,4.80,4.80,5.55,6.10,6.40,0.86,0.86,1.55,2.10,2.10,5.05)
Time_anesthetized.h %>% summary()



colData.df[colData.df$Cell_type=="LHLK"&colData.df$Gender=='Female' &colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000, ] %>% nrow() #41
colData.df[colData.df$Cell_type=="LHLK"&colData.df$Gender=='Female' &colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000, ] %>% rownames() -> g_41LHLK_id
g_41LHLK_colData.df <- colData.df[g_41LHLK_id,]

nrow(colData.df[g_41LHLK_id,])

# State=c(Fed="#DFDF00", Starved="#009F7F" )
g_41LHLK_colData.df %>% select(Cell_id ,Raw_read_num, uniq.reads, Train_protocol, detected.gene.count) -> est.depth
class(est.depth$detected.gene.count)
as.numeric(est.depth$uniq.reads) -> est.depth$uniq.reads
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_scatter_41LHLKs_detected.gene-log.uni.r.08APR2024.pdf",
#     width = 10, height = 8)
sc <- ggplot(est.depth, aes(x = uniq.reads, y = detected.gene.count) )
sc + geom_point(aes(color = Train_protocol), size=5, alpha=.7) +
  scale_color_manual(breaks = c("Fed", "Starved"),values=c("#AAAA00","#00AAAA")) +
  geom_smooth(method = 'lm', se = F, formula = y ~ log(x), size=1.2,color = '#FF7777') +
  # geom_text(aes(label = Cell_id ), size = 6, vjust = 0, hjust = 0, alpha = 0.65,check_overlap = T ) +
  theme_gray(base_size = 20) + scale_x_log10() +
  labs(title=paste("23 Fed LHLKs & 18 Starved LHLKs\nmedian detected gene count ",
                   median(as.numeric(est.depth$detected.gene.count)),"\nmedian unique read count ",
                   median(as.numeric(est.depth$uniq.reads)), sep = ''),
       x ="log10(unique read count)", y = "Detected gene count" ) +
  geom_hline(yintercept= median(as.numeric(est.depth$detected.gene.count)), linetype="dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept= median(as.numeric(est.depth$uniq.reads)), linetype="dashed", color = "black", alpha = 0.5)
dev.off()


## DESeq2 for 18 starved vs 23 fed ####
ncol(count.df)
count.df[,g_41LHLK_id] -> g_41LHLK_count.df
### fwrite the raw count & tpm tab ####
g_41LHLK_count.df[1:17106,] -> g_41LHLK_count.df # remove the rows starting with '__'
# fwrite(data.table(g_41LHLK_count.df,keep.rownames = 'gene_id'),
#           "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/manuscript_G3_2024/Supple/SupplTab2_Patch-seq_41LHLK_count.csv",sep = ',')
# lk_tpm[,.SD,.SDcol=c('gene_id',g_41LHLK_id)] %>% ncol()
# fwrite(lk_tpm[,.SD,.SDcol=c('gene_id',g_41LHLK_id)],
#           "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/manuscript_G3_2024/Supple/SupplTab3_Patch-seq_41LHLK_TPM.csv",sep = ',')


ncol(g_41LHLK_count.df)
g_41LHLK_count.df[rowSums(g_41LHLK_count.df)>0,] %>% nrow()
nrow(g_41LHLK_count.df) 




### detect 12100 genes out of 17106 genes ####

g_41LHLK_colData.df
str(g_41LHLK_colData.df)
class(g_41LHLK_colData.df$Library)
as.factor(g_41LHLK_colData.df[,"Train_protocol"]) -> g_41LHLK_colData.df$Train_protocol
### sanity check heatmap for 23 fed lhlk ####
sanity.c <- c('coGFP','elav','Lk','roX2','repo')
colData.df[c(g_41LHLK_id),c("Library", "Train_protocol")] -> ph.colData
names(ph.colData)[2] <- c("State")
rownames(ph.colData[ph.colData$State=='Fed',]) -> fed_23lhlk_id


lk_tpm.df["Syt1",fed_23lhlk_id]
#       lk1_046 lk1_054 lk1_055 lk1_057 lk1_058 lk3_078 lk3_088 lk3_089
# roX2       0       0       0       0       0       0       0       0
#       lk3_096 lk4_098 lk4_099 lk4_110 lk4_111 lk4_112  lk4_113 lk4_114
# roX2 24.34393       0       0       0       0       0 17.07242       0
#       lk4_115 lk4_116 lk4_117 lk4_118 lk4_119  lk4_121 lk4_122
# roX2       0       0       0       0       0 19.95928       0
lk_tpm.df[c(sanity.c),fed_23lhlk_id] -> sanity.lk
ncol(sanity.lk)

max(log10(sanity.lk+1)) # 4.2
min(log10(sanity.lk+1)) # 0
curBreaks <- seq(0, 4.2, length.out=101)
select(ph.colData, Library) -> ph.colData
ph.clr = list( Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"))

# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap_23fed.lhlk_lib-stat_sanity.log_20240323.pdf", width = 15, height = 3.3)
pheatmap(log10(sanity.lk+1), cluster_rows=F, show_rownames=T,labels_row = c("GFP","elav","Syt1","Lk","roX2","repo"), show_colnames = F, cluster_cols=T, annotation_col= ph.colData, annotation_colors = ph.clr,
         breaks = curBreaks, legend = T, legend_breaks = c(0:3),
         main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 20, treeheight_col = 0, cellwidth = NA)
dev.off()

# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap_23fed.lhlk_stat_sanity.log_07APR2024.pdf", width = 9, height = 3.1)
pheatmap(log10(sanity.lk+1), cluster_rows=F, show_rownames=T,labels_row = c("GFP","elav","Lk","roX2","repo"), show_colnames = F, cluster_cols=T, 
         breaks = curBreaks, legend = T, legend_breaks = c(0:3), legend_labels = c('0','10','100','1000'),
         main = "Transcript abundance", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 20, treeheight_col = 0, cellwidth = NA)
dev.off()
sanity.lk["roX2",] %>% as.numeric() %>% mean()
sanity.lk["repo",] %>% as.numeric() %>% boxplot()

### trsn dive ####
g_41LHLK_id
fed_23lhlk_id
lk_tpm.dt[gene_id=='trsn',.SD,.SDcols = c('gene_id',g_41LHLK_id)]
lk_tpm.dt[gene_id=='trsn',.SD,.SDcols = g_41LHLK_id] %>% t() %>% summary()
lk_tpm.dt[gene_id=='trsn',.SD,.SDcols = fed_23lhlk_id] %>% t() %>% summary()



### Library + Train_protocl only DESeq ####
dds <- DESeqDataSetFromMatrix(countData = g_41LHLK_count.df, 
                              colData = g_41LHLK_colData.df, 
                              design = ~ Library + Train_protocol)
# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
dds.preF
dds["CG32276",]
dds["inc",]
dds["Inos",]
dds.preF["CG32276",]
dds["regucalcin",]
dds["regulocalcin"]
dds.preF["regucalcin",]
dds["CG32276",]
dds.preF["CG32276",]
# Identify the 2000 most variable genes for DE analysis will miss some potential DE genes!! 
assay(dds.preF) %>% log1p %>% rowVars(useNames = T) -> vars
# str(vars)
# names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE) 
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds["Sec61gamma"]
dds.preF["Sec61gamma"]

dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref='Fed')
dds.preF <- DESeq(dds.preF, test="LRT", reduced = ~Library, useT=T, minmu=1e-6,
             minReplicatesForReplace = Inf) # use single-cell recommendation settings
resultsNames(dds.preF)
Lib_res <- results(dds.preF, contrast = c("Train_protocol","Starved","Fed")) # alpha = 0.05, 
Lib_res
summary(Lib_res)
# out of 2000 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 12, 0.6%
# LFC < 0 (down)     : 12, 0.6%
# outliers [1]       : 2, 0.1%
# low counts [2]     : 388, 19%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
Lib_res.Ordered <- Lib_res[order(Lib_res$padj),]
Lib_res.Ordered$log2FoldChange %>% summary()
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -2.326011 -0.229772 -0.004478 -0.009310  0.205097  2.585914 
Lib_res.Ordered$padj <- ifelse(is.na(Lib_res.Ordered$padj), 1, Lib_res.Ordered$padj) # to remove the error message
# Error: logical subscript contains NAs
# write.csv(Lib_res.Ordered, "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/sc.2kVars.DESeq2_41lhlk_lib-state_23MAR2024.csv")
Lib_g41lhlk.gene = rownames(Lib_res.Ordered[Lib_res.Ordered$padj<=0.10,])
Lib_g41lhlk.gene

DESeq2::plotMA(Lib_res, ylim=c(-3,3), ylab="log2 fold change", main='padj<=0.1')
dev.off()
#### heat map ####
colData.df[c(g_41LHLK_id),c("Library", "Train_protocol")] -> ph.colData
names(ph.colData)[2] <- c("State")
rownames(ph.colData[ph.colData$State=='Fed',]) -> fed_23lhlk_id
setdiff(g_41LHLK_id, fed_23lhlk_id) -> starved_18lhlk_id
# brute force to find the col_order
lk_tpm.df[c(Lib_g41lhlk.gene),c(fed_23lhlk_id)] -> de_23.lk
log10(de_23.lk+1)
curBreaks <- seq(0, 3, length.out=101)

ph.clr = list( Gender=c(Female="#FF00FF",Male="#1FDFFF"),
               Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
               State=c(Fed="#DFDF00", Starved="#009F7F" ) )

pheatmap(log10(de_23.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
         main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 7.5, treeheight_row = 0, cellwidth = NA)
dev.off()

lk_tpm.df[c(Lib_g41lhlk.gene),c(starved_18lhlk_id)] -> de_18.lk
log10(de_18.lk+1)
curBreaks <- seq(0, 3, length.out=101)

ph.clr = list( Gender=c(Female="#FF00FF",Male="#1FDFFF"),
               Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
               State=c(Fed="#DFDF00", Starved="#009F7F" ) )

pheatmap(log10(de_18.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
         main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 7.5, treeheight_row = 0, cellwidth = NA)
dev.off()
cluster_order_41lhlk_id = gsub("^","lk",c('4_117','1_046','3_078','4_116','4_098',
                                          '4_099','4_111','3_088','4_113','4_118',
                                          '4_110','4_119','4_114','4_115','4_112',
                                          '4_121','1_054','1_055','1_057','3_096',
                                          '4_122','1_058','3_089',
                                          '3_100','1_043','1_049','1_052','4_028',
                                          '4_105','1_053','1_044','1_050','3_086',
                                          '4_123','3_092','4_108','3_101','4_102',
                                          '4_104','4_103','4_107'
                                          ))
lk_tpm.df[c(Lib_g41lhlk.gene),c(cluster_order_41lhlk_id)] -> de_41.lk
log10(de_41.lk+1)
max(log10(de_41.lk+1)) #3.715
min(log10(de_41.lk+1)) #0
curBreaks <- seq(0, 3.5, length.out=101)
ph.clr = list( #Gender=c(Female="#FF00FF",Male="#1FDFFF"),
               Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
               State=c(Fed="#DFDF00", Starved="#009F7F" ) )

# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap_41lhlk_sc.2kVars.DESeq2-24genes_23MAR2024.pdf", width = 17)
pheatmap(log10(de_41.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T, legend_breaks = c(0:3), 
         main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 15, treeheight_row = 0, cellwidth = NA)
dev.off()
select(ph.colData, State) -> ph.colData
ph.clr[2] -> ph.clr
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap-clusCol_41lhlk_sc.2kVars.DESeq2-24genes-noLib_23MAR2024.pdf", width = 17)
pheatmap(log10(de_41.lk+1), cluster_rows=T, show_rownames=T,show_colnames = F, cluster_cols=T, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T, legend_breaks = c(0:3), legend_labels = c('0','10','100','1000'),
         main = "Transcript abundance", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 15, treeheight_row = 0, cellwidth = NA)
dev.off()
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap_41lhlk_sc.2kVars.DESeq2-24genes-noLib_23MAR2024.pdf", width = 17)
pheatmap(log10(de_41.lk+1), cluster_rows=T, show_rownames=T,show_colnames = F, cluster_cols=F, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T, legend_breaks = c(0:3), legend_labels = c('0','10','100','1000'),
         main = "Transcript abundance", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 15, treeheight_row = 0, cellwidth = NA)
dev.off()
#### relative gene expression ####

# relative gene expression
nrow(de_41.lk)
log2((1 + de_41.lk) / (1 + apply(de_41.lk, 1, mean))) -> rel.de_41.lk


# brute force to find the col_order
rel.de_41.lk[c(Lib_g41lhlk.gene),c(fed_23lhlk_id)] -> de_23.lk
de_23.lk[de_23.lk < -2] <- -2
de_23.lk[de_23.lk > 2] <- 2
curBreaks <- seq(-2, 2, length.out=101)
head(de_23.lk)

pheatmap(de_23.lk, cluster_rows=T, show_rownames=T, cluster_cols=T,  breaks = curBreaks, legend = T,
         main = "Relative Gene Expression (Fed)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=colorRampPalette(c("#998ec3","#f7f7f7","#f1a340"))(100), fontsize = 12, treeheight_row = 0, cellwidth = NA)
dev.off()

rel.de_41.lk[c(Lib_g41lhlk.gene),c(starved_18lhlk_id)] -> de_18.lk
de_18.lk[de_18.lk < -2] <- -2
de_18.lk[de_18.lk > 2] <- 2
curBreaks <- seq(-2, 2, length.out=101)
head(de_18.lk)
pheatmap(de_18.lk, cluster_rows=T, show_rownames=T, cluster_cols=T,  breaks = curBreaks, legend = T,
         main = "Relative Gene Expression (Starved)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=colorRampPalette(c("#998ec3","#f7f7f7","#f1a340"))(100),fontsize = 12, treeheight_row = 0, cellwidth = NA)
dev.off()
cluster_order_41lhlk_id = gsub("^","lk",c('1_054','1_058','4_098','4_099','4_119',
                                          '4_112','4_115','1_055','1_057','1_046',
                                          '3_088','4_121','4_110','3_078','4_114',
                                          '4_116','3_096','4_117','4_118','4_122',
                                          '3_089','4_111','4_113',
                                          '1_043','1_049','1_052','4_105','1_044',
                                          '3_100','4_028','1_050','1_053','3_092',
                                          '3_101','4_108','4_102','4_104','3_086',
                                          '4_123','4_103','4_107'))


log2((1 + de_41.lk) / (1 + apply(de_41.lk, 1, mean))) -> rel.de_41.lk
rel.de_41.lk[rel.de_41.lk < -2] <- -2
rel.de_41.lk[rel.de_41.lk > 2] <- 2
curBreaks <- seq(-2, 2, length.out=101)
head(rel.de_41.lk[,cluster_order_41lhlk_id])

# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.relative_41lhlk_sc.2kVars.DESeq2-24genes-noLib_23MAR2024.pdf")
pheatmap(rel.de_41.lk[,cluster_order_41lhlk_id], cluster_rows=TRUE, show_rownames=T,show_colnames = F,
         cluster_cols=F, annotation_col = ph.colData, annotation_colors = ph.clr, main = "Starvation-up/downregulated genes in LHLK",
         breaks = curBreaks, color = colorRampPalette(c("#998ec3","#f7f7f7","#f1a340"))(100),
         fontsize = 12,# fontsize_col = 20,
         treeheight_row = 0, cellwidth = NA)
dev.off()

## volcano plot ####
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

Lib_res.Ordered

EnhancedVolcano(Lib_res.Ordered,
                lab = rownames(Lib_res.Ordered),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Starved LHLK versus Fed LHLK',
                pCutoff = 1e-1,
                FCcutoff = log2(1.5),
                pointSize = 4.0,
                labSize = 5.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 2/5,
                legendPosition = 'right',
                legendLabSize = 18,
                legendIconSize = 6.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') +
  ggplot2::coord_cartesian(ylim=c(-0.1, 5.8))
dev.off()






# (obsolete) 20 feb vs 18 starved ####
names(lk_tpm.df)[lk_tpm.df["roX2",]==0] -> g_38LHLK_id
g_38LHLK_colData.df <- g_42LHLK_colData.df[g_38LHLK_id,] 

## sanity check heatmap for 20 fed lhlk ####
# to be continued 2024.03.22 ####
sanity.c <- c('GAL4','coGFP','elav', 'Lk', 'trsn','roX2','repo')
lk_tpm.df["roX2",g_38LHLK_id]
lk_tpm.df[c(sanity.c),g_38LHLK_id] -> sanity.lk


ncol(sanity.lk)
log10(sanity.lk+1)
curBreaks <- seq(0, 3, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.lk_sanity.brd.log.20210504.pdf", width = 18, height = 10)
# pheatmap(log10(sanity.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col= NA, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          #          fontsize_col = 20,
#          fontsize = 10, treeheight_row = 0, cellwidth = NA)
# dev.off()




## DESeq2 for 18 starved vs 20 fed ####
ncol(count.df)
count.df[,g_38LHLK_id] -> g_38LHLK_count.df
ncol(g_38LHLK_count.df)
g_38LHLK_colData.df
str(g_38LHLK_colData.df)
class(g_38LHLK_colData.df$Library)
as.factor(g_38LHLK_colData.df[,"Train_protocol"]) -> g_38LHLK_colData.df$Train_protocol
### Library + Train_protocl only DESeq ####
dds <- DESeqDataSetFromMatrix(countData = g_38LHLK_count.df, 
                              colData = g_38LHLK_colData.df, 
                              design = ~ Library + Train_protocol)
keep <- rowSums(counts(dds)) >= 19 # about half of the sample number
dds <- dds[keep,]
dds$Train_protocol <- relevel(dds$Train_protocol, ref='Fed')
dds <- DESeq(dds, test="LRT", reduced = ~Library, useT=T, minmu=1e-6,
             minReplicatesForReplace = Inf) # use single-cell recommendation settings
resultsNames(dds)
Lib_res <- results(dds, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
Lib_res
summary(Lib_res)
Lib_res.Ordered <- Lib_res[order(Lib_res$padj),]
Lib_res.Ordered$log2FoldChange %>% summary()
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -18.83403  -0.21748   0.05247   0.05595   0.30917  18.93562 

Lib_res.Ordered$padj <- ifelse(is.na(Lib_res.Ordered$padj), 1, Lib_res.Ordered$padj) # to remove the error message
# Error: logical subscript contains NAs
# write.csv(Lib_res.Ordered[Lib_res.Ordered$padj<=0.050,], "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/sc.DE_38lhlk_lib-state_23MAR2024.csv")
Lib_g38lhlk.gene = rownames(Lib_res.Ordered[Lib_res.Ordered$padj<=0.050,])
Lib_g38lhlk.gene

DESeq2::plotMA(Lib_res, ylim=c(-6,6), ylab="log2 fold change")
dev.off()
#### heat map ####
colData.df[c(g_38LHLK_id),c("Library", "Train_protocol")] -> ph.colData
names(ph.colData)[2] <- c("State")
nrow(ph.colData)
rownames(ph.colData[ph.colData$State=='Fed',]) -> fed_20lhlk_id
setdiff(g_38LHLK_id, fed_20lhlk_id) -> starved_18lhlk_id
# brute force to find the col_order
lk_tpm.df[c(Lib_g38lhlk.gene),c(fed_20lhlk_id)] -> de_20.lk
log10(de_20.lk+1)
curBreaks <- seq(0, 3, length.out=101)

ph.clr = list( #Gender=c(Female="#FF00FF",Male="#1FDFFF"),
               Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
               State=c(Fed="#DFDF00", Starved="#009F7F" ) )

pheatmap(log10(de_20.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
         main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 7.5, treeheight_row = 0, cellwidth = NA)
dev.off()

# lk_tpm.df[c(Lib_g38lhlk.gene),c(starved_18lhlk_id)] -> de_18.lk
# log10(de_18.lk+1)
# curBreaks <- seq(0, 3, length.out=101)
# 
# ph.clr = list( Gender=c(Female="#FF00FF",Male="#1FDFFF"),
#                Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
#                State=c(Fed="#DFDF00", Starved="#009F7F" ) )
# 
# pheatmap(log10(de_18.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
#          color=viridis_pal(option = "A")(101) ,fontsize = 7.5, treeheight_row = 0, cellwidth = NA)
# dev.off()
cluster_order_38lhlk_id = gsub("^","lk",c(
  '4_119','4_111','4_110','4_122','4_112',
  '4_117','3_078','4_115','4_099','4_118',
  '1_054','3_088','1_046','4_114',
  '4_098','4_116','1_058','3_089','1_055',
  '1_057',
  '1_043','1_052','1_049','3_100','1_053',
  '1_044','1_050','3_086','4_105','4_028',
  '4_123','4_103','4_107','4_104','3_101',
  '4_102','3_092','4_108'))
lk_tpm.df[c(Lib_g38lhlk.gene),c(cluster_order_38lhlk_id)] -> de_38.lk
log10(de_38.lk+1)
max(log10(de_38.lk+1)) #4.068
min(log10(de_38.lk+1)) #0
curBreaks <- seq(0, 3.5, length.out=101)
ph.clr = list( #Gender=c(Female="#FF00FF",Male="#1FDFFF"),
  Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
  State=c(Fed="#DFDF00", Starved="#009F7F" ) )

pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap_38lhlk_de-16genes_23MAR2024.pdf", width = 17)
pheatmap(log10(de_38.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T, legend_breaks = c(0:3), 
         main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 15, treeheight_row = 0, cellwidth = NA)
dev.off()

# lk.4 female only ####
names(lk_tpm.df)[lk_tpm.df["roX2",]==0] -> g_38LHLK_id
g_38LHLK_colData.df <- g_42LHLK_colData.df[g_38LHLK_id,] 

## sanity check heatmap for 20 fed lhlk ####
# to be continued 2024.03.22 ####
sanity.c <- c('GAL4','coGFP','elav', 'Lk', 'trsn','roX2','repo')
lk_tpm.df["roX2",g_38LHLK_id]
lk_tpm.df[c(sanity.c),g_38LHLK_id] -> sanity.lk


ncol(sanity.lk)
log10(sanity.lk+1)
curBreaks <- seq(0, 3, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.lk_sanity.brd.log.20210504.pdf", width = 18, height = 10)
# pheatmap(log10(sanity.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col= NA, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          #          fontsize_col = 20,
#          fontsize = 10, treeheight_row = 0, cellwidth = NA)
# dev.off()




## DESeq2 for 18 starved vs 20 fed ####
ncol(count.df)
count.df[,g_38LHLK_id] -> g_38LHLK_count.df
ncol(g_38LHLK_count.df)
g_38LHLK_colData.df
str(g_38LHLK_colData.df)
class(g_38LHLK_colData.df$Library)
as.factor(g_38LHLK_colData.df[,"Train_protocol"]) -> g_38LHLK_colData.df$Train_protocol
### Library + Train_protocl only DESeq ####
dds <- DESeqDataSetFromMatrix(countData = g_38LHLK_count.df, 
                              colData = g_38LHLK_colData.df, 
                              design = ~ Library + Train_protocol)
keep <- rowSums(counts(dds)) >= 19 # about half of the sample number
dds <- dds[keep,]
dds$Train_protocol <- relevel(dds$Train_protocol, ref='Fed')
dds <- DESeq(dds, test="LRT", reduced = ~Library, useT=T, minmu=1e-6,
             minReplicatesForReplace = Inf) # use single-cell recommendation settings
resultsNames(dds)
Lib_res <- results(dds, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
Lib_res
summary(Lib_res)
Lib_res.Ordered <- Lib_res[order(Lib_res$padj),]
Lib_res.Ordered$log2FoldChange %>% summary()
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -18.83403  -0.21748   0.05247   0.05595   0.30917  18.93562 

Lib_res.Ordered$padj <- ifelse(is.na(Lib_res.Ordered$padj), 1, Lib_res.Ordered$padj) # to remove the error message
# Error: logical subscript contains NAs
# write.csv(Lib_res.Ordered[Lib_res.Ordered$padj<=0.050,], "/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/sc.DE_38lhlk_lib-state_23MAR2024.csv")
Lib_g38lhlk.gene = rownames(Lib_res.Ordered[Lib_res.Ordered$padj<=0.050,])
Lib_g38lhlk.gene

DESeq2::plotMA(Lib_res, ylim=c(-6,6), ylab="log2 fold change")
dev.off()
#### heat map ####
colData.df[c(g_38LHLK_id),c("Library", "Train_protocol")] -> ph.colData
names(ph.colData)[2] <- c("State")
nrow(ph.colData)
rownames(ph.colData[ph.colData$State=='Fed',]) -> fed_20lhlk_id
setdiff(g_38LHLK_id, fed_20lhlk_id) -> starved_18lhlk_id
# brute force to find the col_order
lk_tpm.df[c(Lib_g38lhlk.gene),c(fed_20lhlk_id)] -> de_20.lk
log10(de_20.lk+1)
curBreaks <- seq(0, 3, length.out=101)

ph.clr = list( #Gender=c(Female="#FF00FF",Male="#1FDFFF"),
  Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
  State=c(Fed="#DFDF00", Starved="#009F7F" ) )

pheatmap(log10(de_20.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
         main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 7.5, treeheight_row = 0, cellwidth = NA)
dev.off()

# lk_tpm.df[c(Lib_g38lhlk.gene),c(starved_18lhlk_id)] -> de_18.lk
# log10(de_18.lk+1)
# curBreaks <- seq(0, 3, length.out=101)
# 
# ph.clr = list( Gender=c(Female="#FF00FF",Male="#1FDFFF"),
#                Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
#                State=c(Fed="#DFDF00", Starved="#009F7F" ) )
# 
# pheatmap(log10(de_18.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
#          color=viridis_pal(option = "A")(101) ,fontsize = 7.5, treeheight_row = 0, cellwidth = NA)
# dev.off()
cluster_order_38lhlk_id = gsub("^","lk",c(
  '4_119','4_111','4_110','4_122','4_112',
  '4_117','3_078','4_115','4_099','4_118',
  '1_054','3_088','1_046','4_114',
  '4_098','4_116','1_058','3_089','1_055',
  '1_057',
  '1_043','1_052','1_049','3_100','1_053',
  '1_044','1_050','3_086','4_105','4_028',
  '4_123','4_103','4_107','4_104','3_101',
  '4_102','3_092','4_108'))
lk_tpm.df[c(Lib_g38lhlk.gene),c(cluster_order_38lhlk_id)] -> de_38.lk
log10(de_38.lk+1)
max(log10(de_38.lk+1)) #4.068
min(log10(de_38.lk+1)) #0
curBreaks <- seq(0, 3.5, length.out=101)
ph.clr = list( #Gender=c(Female="#FF00FF",Male="#1FDFFF"),
  Library=c(Lk.1="#fbb4ae", Lk.3="#b3cde3", Lk.4="#ccebc5"),
  State=c(Fed="#DFDF00", Starved="#009F7F" ) )

pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap_38lhlk_de-16genes_23MAR2024.pdf", width = 17)
pheatmap(log10(de_38.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=F, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T, legend_breaks = c(0:3), 
         main = "Marker Gene Expression (log10 TPM)", display_numbers = F, number_color = "grey10",angle_col = 90,
         color=viridis_pal(option = "A")(101) ,fontsize = 15, treeheight_row = 0, cellwidth = NA)
dev.off()

# to be continued 2024.03.22 ####








# # If used in published research, please cite:
# #   Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
# # sequence count data: removing the noise and preserving large differences.
# # Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
# resLFC <- lfcShrink(dds, coef = "Train_protocol_Starved_vs_Fed",
#                     type="apeglm")
# resLFC
# plotMA(resLFC, ylim=c(-2,2))












# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
filter <- rowSums(assay(dds)>1)>5
table(filter) # FALSE 10276, TRUE 6830
# Identify the 2000 most variable genes
assay(dds.preF) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE)
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
dds.preF.analysis <- DESeq(dds.preF)
res5.S2F <- results(dds.preF.analysis, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
summary(res5.S2F)
# 1 up-regulated and 3 down-regulated DE genes (padj < 0.05)
res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
res5.S2F.Ordered
# write.csv(as.data.frame(res5.S2F.Ordered),
#           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk_L.9Starved-13Fed_DESeq2_20201020.csv")
res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])
S2F.gene
data.table(gene = rownames(as.data.frame(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])),
           as.data.frame(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,]),
           method = 'DESeq2', dataset = 'Lk.1,Lk.3,Lk.4', stage = 'early' ) -> L.all.DESeq2.de.tab
L.all.edgeR.de.tab







## 24 Fed vs 18 Starved ####
# ERCCoverTranscript_ratio<0.09 &
# Spike.in.linearity>0.80 &
# Assigned.reads >250000
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000 & 
             colData.df$Cell_id!=45, "detected.gene.count"] %>% mean() # 5985.262
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000 & 
             colData.df$Cell_id!=45, "detected.gene.count"] %>% median() # 5821
rownames(g.colData.df) -> g.lk.id
lk1.counts[lk3.counts, on=c("gene_id==gene_id")] -> lk_counts
lk_counts[lk4.counts, on=c("gene_id==gene_id")] -> lk_counts
f_dowle2(lk_counts)
data.frame(lk_counts[,2:(length(c(lk1_id, lk3_id, lk4_id))+1)], row.names = lk_counts$gene_id) -> count.df

count.df[,g.lk.id] -> g.count.df
ncol(g.count.df)


## Heatmap for publication ####
sanity.c <- c('GAL4','coGFP','roX2','elav','repo','Lk','trsn','Gad1')
# all cells
lk_tpm.df[c(sanity.c),] -> sanity.lk
log10(sanity.lk+1)
curBreaks <- seq(0, 3, length.out=101)
colData.df[,c("Cell_type", "Train_protocol")] -> ph.colData
ph.colData["lk1_041",1]<- c("MBON.g3bp1")
names(ph.colData) <- c('Cell_type',"State")
ph.clr = list( Cell_type=c(MBON.g3bp1="#FF00FF",LHLK="#1FDFFF"), State=c(Fed="#DFDF00", Starved="#009F7F" ) )

# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.lk.all-LHLKs_sanity.log.20210504.pdf", width = 17)
# pheatmap(log10(sanity.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col = ph.colData, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          fontsize = 7.5, treeheight_row = 0, cellwidth = NA)
# dev.off()

# 42 LHLKs
lk_tpm.df[c(sanity.c),g.lk.id] -> sanity.lk
log10(sanity.lk+1)
curBreaks <- seq(0, 3, length.out=101)
data.frame(State= colData.df[g.lk.id, "Train_protocol"], row.names=g.lk.id ) -> ph.colData.0
ph.clr = list( State=c(Fed="#DFDF00", Starved="#009F7F") )
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.lk.42LHLKs_sanity.log.20210504.pdf", width = 14)
# pheatmap(log10(sanity.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col = ph.colData.0, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          fontsize = 8.5, treeheight_row = 0, cellwidth = NA)
# dev.off()


## Naive transcriptome ####
intersect(rownames(colData.df[colData.df$Train_protocol=='Fed',]), g.lk.id) -> naive.id
length(naive.id) # 24 naive cells
data.frame(State= colData.df[naive.id, "Train_protocol"], row.names=naive.id) -> ph.colData
str(ph.colData)
ph.clr = list( State=c(Fed="#DFDF00", Starved="#009F7F") )
gene.list <- c('Lk','trsn','ap','AMPKalpha')
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
log10(naive.lk+1)
curBreaks <- seq(0, 3, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.24Fed.lk_LHLK.markers.log.20210504.pdf", height = 4, width = 12)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()


# Neurotransmitter
fread("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/NGS/annotation_gtf/Gene.list.Mx20200701.tsv") -> GOI
GOI[Group=="Neurotransmitter",Symbol] -> gene.list
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
log10(naive.lk+1)
curBreaks <- seq(0, 2.5, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.24Fed.lk_nt.log.20210504.pdf", width = 12)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=T, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Neurotransmitter Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()

# Neurotransmitter receptor
GOI[Group=="Neurotransmitter receptor",Symbol] -> gene.list
lk_tpm.df[rownames(lk_tpm.df)%in%gene.list ,naive.id] -> naive.lk
# missing "Grik" and "Octalpha2R" due to the outdated annotation file genereated in 2017.
log10(naive.lk+1)
curBreaks <- seq(0, 3, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.24Fed.lk_ntR.log.20210504.pdf", height = 13, width = 12)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=T, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Neurotransmitter Receptor Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()

# Heatmap for NP
np <- c("Akh","amn","AstA","AstC","AstCC","Burs","Capa","CCAP","CCHa1","CCHa2","CNMa","Crz","Dh31","Dh44","Dsk","Eh",
        "ETH","FMRFa","Gpa2","Gpb5","Hug","Ilp1","Ilp2","Ilp3","Ilp4","Ilp5","Ilp6","Ilp7","Ilp8","ITP","Lk",
        "Mip","Ms","NPF","Nplp1","Nplp3","Nplp4","Orcokinin","Pburs","Pdf","Proc","Ptth","RYa","SIFa",
        "sNPF","SP","spab","Tk")

gene.list <- c(np)
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
lk_tpm.df["AstCC",]
log10(naive.lk+1)
curBreaks <- seq(0, 3, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.24Fed.lk_np.log.20210504.pdf", height = 12, width = 12)
# pheatmap(log10(naive.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Neuropeptide Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()
# NP receptor
GOI[150:200,]
GOI[Group=="Neuropeptide receptor",Symbol] -> gene.list
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
log10(naive.lk+1)
curBreaks <- seq(0, 2.4, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.24Fed.lk_NpR.log.20210504.pdf", height = 12, width = 12)
# pheatmap(log10(naive.lk+1), cluster_rows=T, show_rownames=T, cluster_cols=T, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Neuropeptide Receptor Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()

# Gap junctions
GOI[Group=="gap junction",Symbol] -> gene.list
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
log10(naive.lk+1)
curBreaks <- seq(0, 2.5, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.24Fed.lk_GJ.log.20210504.pdf", width = 12)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Gap Junction Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()


# Size estimation & power analysis ####
## 20201118 # to be continued
# Use only lk.4 as prior dataset!
g.colData.df[c(g.lk.id),"Train_protocol"] -> lk.index
gsub("Starved",1,lk.index) -> lk.index
gsub("Fed",0,lk.index) -> lk.index
class(lk.index)
dataMatrix = g.count.df
#Estimate the gene read count and dispersion distribution
dataMatrixDistribution<-est_count_dispersion(dataMatrix, group=c(as.numeric(lk.index)), subSampleNum = length(lk.index))
# Disp = 0.23399 , BCV = 0.4837 # "2020-11-18 14:03:26 EST" from lk.1,3,4
# Disp = 0.22491 , BCV = 0.4743 # "2021-05-04 13:36:04 EDT" from lk1,3,4
# Get timestamp from system time
# Sys.time()
selectedGenes = c("fit","inc")
g.count.df[selectedGenes,]
# Test  Number of samples in each group: 17;
#       Minimal fold change between two groups: 2;
#       False discovery rate: 0.01;
#       only 'fit' and 'inc' two genes
powerDistribution<-est_power_distribution(n=34,f=0.01,rho=2, # 17*2 = 34 is the number of valid primers
                                          distributionObject=dataMatrixDistribution,
                                          selectedGenes=selectedGenes,
                                          storeProcess=TRUE)
str(powerDistribution)
mean(powerDistribution$power)
# [1] 0.2426721 from 17 samples in each group
# 0.7436811 from 34 samples in each group
# 0.8425432 probability to find the significant genes, fit and inc, from 40 samples in each group


# Estimate sample size 
#   with 0.8 power, 0.01 false discovery rat, lk.1,3,4 as prior dataset

sample_size_distribution(power=0.8,f=0.01,distributionObject=dataMatrixDistribution,
                         selectedGenes=selectedGenes,
                         repNumber=100,
                         showMessage=TRUE)
# 38 samples in each group for 'fit' and 'inc'
# 25 samples in each group for random genes


dataMatrix[,grep("lk4", g.lk.id)] -> dataMatrix.lk4
lk.index[grep("lk4", g.lk.id)] -> lk.index.lk4
dataMatrixDistribution_lk4<-est_count_dispersion(dataMatrix.lk4, group=c(as.numeric(lk.index.lk4)), subSampleNum = length(lk.index.lk4))
# Disp = 0.19475 , BCV = 0.4413  # "2021-05-04 13:39:42 EDT" from lk.4
# Test  Number of samples in each group: 17;
#       Minimal fold change between two groups: 2;
#       False discovery rate: 0.01;
#       only 'fit' and 'inc' two genes
powerDistribution<-est_power_distribution(n=17,f=0.01,rho=2, # 17*2 = 34 is the number of valid primers
                                          distributionObject=dataMatrixDistribution_lk4,
                                          selectedGenes=selectedGenes,
                                          storeProcess=TRUE)
str(powerDistribution)
mean(powerDistribution$power)
# 0.3470254 probability to find the significant genes, fit and inc, from 17 samples in each group
# 0.8627373 probability to find the significant genes, fit and inc, from 34 samples in each group
sample_size_distribution(power=0.8,f=0.01,distributionObject=dataMatrixDistribution_lk4,
                         # selectedGenes=selectedGenes,
                         # repNumber=100,
                         showMessage=TRUE)
# 32 samples in each group for 'fit' and 'inc'
# 27 samples in each group for random genes



# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_bar_42.lk_de.g.20210504.pdf", width = 12)
# ggplot(data = g.colData.df, aes(x = rownames(g.colData.df), y=detected.gene.count) ) + geom_bar(stat = "identity", width = 0.6) +
#   scale_x_discrete(limits=rownames(g.colData.df) ) + scale_y_continuous(labels=c("0","2k","4k","6k","8k") ) + theme_gray(base_size = 16) +
#   theme(axis.text.x = element_text(angle=45 , hjust = 1, vjust = 1) ) +
#   labs(title="24 Fed vs 18 Starved LHLKs", x ="Cell_id", y = "Detected gene count") + geom_hline(yintercept=median(as.vector(colSums(g.count.df>0))), linetype="dashed", color = "blue") +
#   geom_text(x = 36.5, y = 7600,label = paste('median = ',median(as.vector(colSums(g.count.df>0))), sep = ''), color="blue", size = 8, family="Courier", fontface="plain")
# dev.off()


# DE analysis ####
## edgeR ####
## Test whether CEL-Seq data is zero-inflated using Scatterplots of 
##   the estimated biological coefficient of variation (BCV, defined 
##   as the square root of the negative binomial dispersion parameter φ)
##   against average log counts per million (CPM) computed using edgeR
nrow(g.colData.df)
ncol(g.count.df)

c(as.vector(g.colData.df$Train_protocol) ) -> group
y <- DGEList(counts=g.count.df, group=group)
y$samples
y$counts['coGFP',]

# Filtering
keep <- rowSums(cpm(y)>1) >= 21 # Total 42 samples, use 21 as cutoff
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)
# 6103 genes passed

# Normalization for RNA composition
y <- calcNormFactors(y)
y$samples

# group lib.size norm.factors
# lk1_043 Starved    37888    1.0817666
# lk1_044 Starved    63877    0.9949394
# lk1_046     Fed    65016    0.9536273
# lk1_049 Starved    32038    0.9937821
# lk1_050 Starved    51291    0.9234947
# lk1_052 Starved    60141    1.1427772
# lk1_053 Starved    67708    0.8891842
# lk1_054     Fed    45083    0.9024079
# lk1_055     Fed    24212    1.0742303
# lk1_057     Fed    58077    1.0080940
# lk1_058     Fed    77167    0.9585607
# lk3_078     Fed    33034    0.9426423
# lk3_086 Starved    46519    1.1035288
# lk3_088     Fed    49514    0.9027646
# lk3_089     Fed    34597    1.0476552
# lk3_092 Starved    59485    0.9197285
# lk3_096     Fed    38764    1.0752256
# lk3_100 Starved    37226    1.0666210
# lk3_101 Starved    47194    0.9521675
# lk4_028 Starved    84202    1.0626977
# lk4_029     Fed    31296    1.0874000
# lk4_098     Fed    32335    0.9302239
# lk4_099     Fed   117425    1.0563565
# lk4_102 Starved    42361    0.9559104
# lk4_103 Starved   118060    1.0667691
# lk4_104 Starved    45748    0.9446106
# lk4_105 Starved    37125    1.0144414
# lk4_107 Starved   144237    1.0354262
# lk4_108 Starved    74884    0.9477132
# lk4_110     Fed    65061    0.9675223
# lk4_111     Fed    99703    1.0491896
# lk4_112     Fed    77678    0.9383624
# lk4_113     Fed    55402    1.0818435
# lk4_114     Fed    70727    0.9140284
# lk4_115     Fed    34970    0.9630309
# lk4_116     Fed    36334    0.9222424
# lk4_117     Fed    35931    1.0124447
# lk4_118     Fed    37418    1.0868213
# lk4_119     Fed    41198    1.0061766
# lk4_121     Fed    47149    0.9961685
# lk4_122     Fed    51053    1.0867592
# lk4_123 Starved    92599    1.0338554

# plotBCV after estimating common dispersion and tagwise dispersions
#  in one run (recommended):
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)
pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/BCV.plot_42gLHLK_Train.Prot_20210504.pdf")
plotBCV(y)
dev.off()

et <- exactTest(y, pair = c("Fed","Starved"))
topTags(et, n = 20)
summary(decideTests(et)) # at 5% FDR
# 3 down & 8 up genes coming out of default edgeR analysis from Starved - Fed comparison
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
topTags(lrt, n = 35)
# Coefficient:  groupStarved 
# logFC    logCPM       LR       PValue          FDR
# fit       -2.0793346  8.143579 37.15588 1.090533e-09 6.655526e-06
# CG5773    -1.8575740  7.632175 32.07765 1.481318e-08 3.728368e-05
# Got2       1.5236413  9.005843 31.66417 1.832722e-08 3.728368e-05
# Yp2       -1.2040163 12.489523 29.75807 4.894646e-08 7.468007e-05
# AttA       2.7322735  7.006899 21.10072 4.357648e-06 5.318946e-03
# Cyp28d1    1.2332024  7.589718 17.91108 2.314692e-05 2.211949e-02
# Ets97D    -2.1243513  5.605779 17.73655 2.537054e-05 2.211949e-02
# AGBE       1.2536807  6.682111 17.02509 3.688923e-05 2.814187e-02
# CG6767     0.7495242  9.797972 16.67417 4.438110e-05 3.009532e-02
# AttC       2.4659835  7.208748 16.45995 4.968868e-05 3.032500e-02
# Scsalpha   0.8735460  8.835456 15.95279 6.494211e-05 3.603107e-02
# rho        1.3378598  6.798248 15.69753 7.432109e-05 3.779847e-02
# CR44396   -1.8195623  5.799212 15.50588 8.224881e-05 3.861265e-02
# CG3348    -2.5989777  6.032461 14.85005 1.164052e-04 4.820809e-02
# AdSL       1.3520049  6.667312 14.77658 1.210294e-04 4.820809e-02
# CG31075    1.3974071  6.672788 14.69494 1.263853e-04 4.820809e-02
# inc       -1.0923609  7.206975 14.35819 1.511208e-04 5.425237e-02
# Sodh-1     1.5784045  8.134049 14.08256 1.749585e-04 5.932065e-02
rownames(topTags(lrt, n = 35,p.value = 0.05)) -> edgeR.DE.gene
summary(decideTests(lrt))
# 6 down-regulated genes & 10 up-regulated genes coming out of 
#the likelihood ratio test edgeR analysis from Starved - Fed comparison
# write.csv(as.data.frame(topTags(lrt,n = 6103)),
#           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/lk_18Starved-24Fed_edgeR_20210504.csv")
# To fill in the missing value # 20201118
data.table(gene = rownames(as.data.frame(topTags(lrt,n = 6103))),as.data.frame(topTags(lrt,n = 6103)),
           method = 'edgeR.lrt', dataset = 'Lk1,Lk3,Lk4' ) -> all.gene_edgeR.24f18s.de.dt

data.table(gene = rownames(as.data.frame(topTags(lrt,n = 16))),as.data.frame(topTags(lrt,n = 16)),
           method = 'edgeR.lrt', dataset = 'Lk1,Lk3,Lk4' ) -> all.edgeR.de.tab

## DESeq2 ####
g.colData.df
str(g.colData.df)
as.factor(g.colData.df[,"Train_protocol"]) -> g.colData.df$Train_protocol
dds <- DESeqDataSetFromMatrix(countData = g.count.df, 
                              colData = g.colData.df, 
                              design = ~ Train_protocol)
# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
filter <- rowSums(assay(dds)>1)>5
table(filter) # FALSE 10281, TRUE 6830
# Identify the 2000 most variable genes for DE analysis will miss some potential DE genes!! 
assay(dds.preF) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE)
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
dds.preF.analysis <- DESeq(dds.preF)
res5.S2F <- results(dds.preF.analysis, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
summary(res5.S2F)
# 5 up-regulated and 5 down-regulated DE genes (padj < 0.05)
res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
res5.S2F.Ordered
write.csv(as.data.frame(res5.S2F.Ordered),
          file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/lk_18Starved-24Fed_DESeq2_20210504.csv")
res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])
S2F.gene


data.table(gene = rownames(as.data.frame(res5.S2F.Ordered)),
           as.data.frame(res5.S2F.Ordered),
           method = 'DESeq2', dataset = 'Lk1,Lk3,Lk4' ) -> all.gene_DESeq2.24f18s.de.dt
all.gene_DESeq2.24f18s.de.dt[padj<0.05,]
all.gene_DESeq2.24f18s.de.dt
all.gene_edgeR.24f18s.de.dt
merge(all.gene_DESeq2.24f18s.de.dt, all.gene_edgeR.24f18s.de.dt, by = c('gene'), all = T, suffixes = c("_DESeq2", "_edgeR") ) ->all.gene_24f18s.de.dt 


data.table(gene = rownames(as.data.frame(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])),
           as.data.frame(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,]),
           method = 'DESeq2', dataset = 'Lk1,Lk3,Lk4' ) -> all.DESeq2.de.tab
all.edgeR.de.tab


## Make a bar chart for the DE genes from edgeR & DESeq2 ####
g.count.df
g.colData.df
bar.order = c(rownames(g.colData.df[g.colData.df$Train_protocol=='Fed',]), rownames(g.colData.df[g.colData.df$Train_protocol=='Starved',]) )
edgeR.DE.gene
unique(c(edgeR.DE.gene,S2F.gene)) -> DE.gene
lk_tpm.df[rownames(lk_tpm.df)%in%DE.gene , rownames(g.colData.df)] -> gene.matrix
lk_tpm.df['CG33926',]
str(gene.matrix)
length(DE.gene)
nrow(gene.matrix)
mx = max(gene.matrix)*4/5
for (i in 1:nrow(gene.matrix)){
  gather(gene.matrix[DE.gene[i],]) -> df
  # mx = max(df)*4/5
  pdf(paste("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/barplot.",DE.gene[i],".42lhlk.20210504.pdf", sep = "") )
  print(ggplot(df, aes(x = key, y=value)) + geom_bar(stat = "identity", width = 0.75) +
          theme_gray(base_size = 10) + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1 )) +
          labs(title = "24 Fed vs 18 Starved",x ="sample", y = paste("tpm of ",DE.gene[i], sep = "")) + geom_vline(xintercept=c(24.5), linetype="dashed", color = "red") +
          scale_x_discrete(limits=bar.order) #+
       #geom_text(x = 3, y = mx ,label = "Fed", color="light brown", size = 9, family="Courier", fontface="plain")
        # geom_text(x = 11, y = max(df)*3/4 ,label = "Starved", color="dark green", size = 9, family="Courier", fontface="plain")
  )
  dev.off()
}



## Make a bar chart for the DE genes exclusively from DESeq2 ####
g.count.df
g.colData.df
bar.order = c(rownames(g.colData.df[g.colData.df$Train_protocol=='Fed',]), rownames(g.colData.df[g.colData.df$Train_protocol=='Starved',]) )
lk_tpm.df[rownames(lk_tpm.df)%in%S2F.gene , rownames(g.colData.df)] -> gene.matrix
lk_tpm.df['CG33926',]
str(gene.matrix)
length(S2F.gene)
nrow(gene.matrix)
mx = max(gene.matrix)*4/5
for (i in 1:nrow(gene.matrix)){
  gather(gene.matrix[S2F.gene[i],]) -> df
  # mx = max(df)*4/5
  pdf(paste("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/barplot.DESeq2_",S2F.gene[i],".42lhlk.20210504.pdf", sep = "") )
  print(ggplot(df, aes(x = key, y=value)) + geom_bar(stat = "identity", width = 0.75) + 
          theme_gray(base_size = 10) + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1 )) +
          labs(title = paste("DESeq2 analyzing 24 Fed vs 18 Starved -- ",S2F.gene[i],sep = "" ),x ="sample", y = paste("tpm of ",S2F.gene[i], sep = "")) + geom_vline(xintercept=c(24.5), linetype="dashed", color = "red") +
          scale_x_discrete(limits=bar.order) #+
  )
  dev.off()
}




# Lk.4 only ####
# colData.df$Library
lk4_id
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000 & 
             colData.df$Cell_id!=45 & colData.df$Library =="Lk.4", ] -> g.colData.df # In DPM project, 0.035 and 0.8 as thresholds
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000 & 
             colData.df$Cell_id!=45 & colData.df$Library =="Lk.4", "Train_protocol"] %>% table()
# 15 Fed vs 8 Starved 
rownames(setdiff(colData.df,g.colData.df)) # Remove "lk.4_025" "lk.4_059" "lk.4_106" "lk.4_109" "lk.4_120"

colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000 &
             colData.df$Cell_id!=45 & colData.df$Library =="Lk.4", "detected.gene.count"] %>% mean() # 6121.522
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000 &
             colData.df$Cell_id!=45 & colData.df$Library =="Lk.4", "detected.gene.count"] %>% median() # 5950
rownames(g.colData.df) -> g.lk.id
data.frame(lk_counts[,2:(length(c(lk1_id, lk3_id, lk4_id))+1)], row.names = lk_counts$gene_id) -> count.df
count.df[,g.lk.id] -> g.count.df
ncol(g.count.df)

# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/plot_bar_23.lk4_de.g.20210504.pdf")
# ggplot(data = g.colData.df, aes(x = rownames(g.colData.df), y=detected.gene.count) ) + geom_bar(stat = "identity", width = 0.6) +
#   scale_x_discrete(limits=rownames(g.colData.df) ) + scale_y_continuous(labels=c("0","2k","4k","6k","8k") ) + theme_gray(base_size = 16) +
#   theme(axis.text.x = element_text(angle=45 , hjust = 1, vjust = 1) ) +
#   labs(title="15 Fed vs 8 Starved LHLKs in Lk.4", x ="Cell_id", y = "Detected gene count") + geom_hline(yintercept=median(as.vector(colSums(g.count.df>0))), linetype="dashed", color = "blue") +
#   geom_text(x = 19, y = 7000,label = paste('median',median(as.vector(colSums(g.count.df>0))), sep = ' '), color="blue", size = 6, family="Courier", fontface="plain")
# dev.off()

# Heatmap for publication
sanity.c <- c('GAL4','coGFP','roX2','elav','repo','Lk','trsn','Gad1')

# 23 LHLKs from Lk.4
lk_tpm.df[c(sanity.c),g.lk.id] -> sanity.lk
log10(sanity.lk+1)
curBreaks <- seq(0, 3, length.out=101)
data.frame(State= colData.df[g.lk.id, "Train_protocol"], row.names=g.lk.id) -> ph.colData.0
ph.clr = list( State=c(Fed="#DFDF00", Starved="#009F7F") )
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/heatmap.lk4.23LHLKs_sanity.log.20210504.pdf", height = 4, width = 10)
# pheatmap(log10(sanity.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col = ph.colData.0, annotation_colors = ph.clr, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",angle_col = 90,
#          fontsize = 8.5, treeheight_row = 0, cellwidth = NA)
# dev.off()


## Naive transcriptome
intersect(rownames(colData.df[colData.df$Train_protocol=='Fed',]), g.lk.id) -> naive.id
length(naive.id) # 15 naive cells
data.frame(State= colData.df[naive.id, "Train_protocol"], row.names=naive.id) -> ph.colData
str(ph.colData)
ph.clr = list( State=c(Fed="#DFDF00", Starved="#009F7F") )
gene.list <- c('Lk','trsn','ap','AMPKalpha')
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
log10(naive.lk+1)
curBreaks <- seq(0, 3, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/heatmap.15Fed.lk4_LHLK.markers.log.20201013.pdf", height = 4, width = 10)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Marker Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()


# Neurotransmitter
fread("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/NGS/annotation_gtf/Gene.list.Mx20200701.tsv") -> GOI
GOI[Group=="Neurotransmitter",Symbol] -> gene.list
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
log10(naive.lk+1)
curBreaks <- seq(0, 2.5, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/heatmap.15Fed.lk4_nt.log.20201013.pdf", width = 10)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=T, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Neurotransmitter Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()

# Neurotransmitter receptor
GOI[Group=="Neurotransmitter receptor",Symbol] -> gene.list
lk_tpm.df[rownames(lk_tpm.df)%in%gene.list ,naive.id] -> naive.lk
# missing "Grik" and "Octalpha2R" due to the outdated annotation file genereated in 2017.
log10(naive.lk+1)
curBreaks <- seq(0, 3, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/heatmap.15Fed.lk4_ntR.log.20201013.pdf", height = 13, width = 10)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=T, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Neurotransmitter Receptor Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()

# Heatmap for NP
np <- c("Akh","amn","AstA","AstC","AstCC","Burs","Capa","CCAP","CCHa1","CCHa2","CNMa","Crz","Dh31","Dh44","Dsk","Eh",
        "ETH","FMRFa","Gpa2","Gpb5","Hug","Ilp1","Ilp2","Ilp3","Ilp4","Ilp5","Ilp6","Ilp7","Ilp8","ITP","Lk",
        "Mip","Ms","NPF","Nplp1","Nplp3","Nplp4","Orcokinin","Pburs","Pdf","Proc","Ptth","RYa","SIFa",
        "sNPF","SP","spab","Tk")

gene.list <- c(np)
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
lk_tpm.df["AstCC",]
log10(naive.lk+1)
curBreaks <- seq(0, 3, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/heatmap.15Fed.lk4_np.log.20201013.pdf", height = 12, width = 10)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=T, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Neuropeptide Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()
# NP receptor
GOI[150:200,]
GOI[Group=="Neuropeptide receptor",Symbol] -> gene.list
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
log10(naive.lk+1)
curBreaks <- seq(0, 2.4, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/heatmap.15Fed.lk4_NpR.log.20201013.pdf", height = 12, width = 10)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=T, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Neuropeptide Receptor Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()

# Gap junctions
GOI[Group=="gap junction",Symbol] -> gene.list
lk_tpm.df[c(gene.list),naive.id] -> naive.lk
log10(naive.lk+1)
curBreaks <- seq(0, 2.5, length.out=101)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/heatmap.15Fed.lk4_GJ.log.20201013.pdf", width = 10)
# pheatmap(log10(naive.lk+1), cluster_rows=F, show_rownames=T, cluster_cols=F, annotation_col = NA, annotation_colors = NA, breaks = curBreaks, legend = T,
#          main = "Gap Junction Gene Expression (log10 TPM)", display_numbers = T, number_color = "grey10",
#          fontsize = 11, treeheight_row = 0, cellwidth = NA)
# dev.off()

## DE analysis
# edgeR
## Test whether CEL-Seq data is zero-inflated using Scatterplots of 
##   the estimated biological coefficient of variation (BCV, defined 
##   as the square root of the negative binomial dispersion parameter φ)
##   against average log counts per million (CPM) computed using edgeR
nrow(g.colData.df)
ncol(g.count.df)

c(as.vector(g.colData.df$Train_protocol) ) -> group
y <- DGEList(counts=g.count.df, group=group)
y$samples
y$counts['coGFP',]

# Filtering
keep <- rowSums(cpm(y)>1) >= 12 # Total 23 samples, use 12 as cutoff
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)
# 6211 genes passed

# Normalization for RNA composition
y <- calcNormFactors(y)
y$samples

# group lib.size norm.factors
# lk4_028 Starved    84427    1.0223547
# lk4_029     Fed    31661    1.0852961
# lk4_098     Fed    32410    0.9486539
# lk4_099     Fed   117632    0.9914865
# lk4_102 Starved    42451    1.0056718
# lk4_103 Starved   118337    1.0158238
# lk4_104 Starved    45856    1.0099677
# lk4_105 Starved    37191    0.9868167
# lk4_107 Starved   144541    1.0489324
# lk4_108 Starved    75058    0.9240976
# lk4_110     Fed    65177    0.9769846
# lk4_111     Fed    99951    1.0250048
# lk4_112     Fed    77820    0.9359905
# lk4_113     Fed    55511    0.9857444
# lk4_114     Fed    70861    0.9534364
# lk4_115     Fed    35077    1.0133587
# lk4_116     Fed    36393    0.9463267
# lk4_117     Fed    36016    1.0153783
# lk4_118     Fed    37508    1.0841010
# lk4_119     Fed    41301    1.0264594
# lk4_121     Fed    47258    0.9954498
# lk4_122     Fed    51156    0.9805681
# lk4_123 Starved    92921    1.0420526

# plotBCV after estimating common dispersion and tagwise dispersions
#  in one run (recommended):
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/BCV.plot_23gLHLK-Lk4_Train.Prot_20201013.pdf")
# plotBCV(y)
# dev.off()

et <- exactTest(y, pair = c("Fed","Starved"))
topTags(et, n = 20)
summary(decideTests(et)) # at 5% FDR
# 5 down & 20 up genes coming out of default edgeR analysis from Starved - Fed comparison
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
topTags(lrt, n = 35)
# Coefficient:  groupStarved 
# logFC    logCPM       LR       PValue          FDR
# AttC        2.8186444  7.551535 49.97625 1.556179e-12 9.665429e-09
# AttB        2.6623026  8.490016 39.68705 2.980953e-10 7.214132e-07
# AttA        3.2559775  7.436179 39.38224 3.484527e-10 7.214132e-07
# Sodh-1      2.2672323  7.779799 38.06477 6.843472e-10 1.062620e-06
# CG16772     2.0036435  8.298828 34.68178 3.882486e-09 4.822824e-06
# Got2        1.6705493  8.989480 33.25109 8.099423e-09 8.384253e-06
# fit        -2.0977029  8.049044 27.11529 1.916763e-07 1.700717e-04
# PGRP-SB1    2.1368678  7.123530 26.06543 3.300413e-07 2.562358e-04
# Spn28Dc     1.7357653  8.357051 24.63659 6.922495e-07 4.341112e-04
# DptB        3.3953112  6.898850 24.61805 6.989392e-07 4.341112e-04
# GNBP-like3  2.4883146  7.282809 20.77853 5.155761e-06 2.911130e-03
# Nmdmc       1.6832277  7.066110 20.48369 6.014140e-06 3.112818e-03
# CG5773     -1.9426268  7.549966 20.15992 7.122972e-06 3.403137e-03
# CG13607    -2.5132736  6.632399 18.77666 1.469544e-05 6.519525e-03
# Yp2        -1.3109672 12.550220 18.41080 1.780456e-05 7.372273e-03
# Cyp28d1     1.4126921  7.408610 17.84443 2.397193e-05 9.305604e-03
# CecA1       2.1981619  7.653763 16.25787 5.527953e-05 2.019654e-02
# CecA2       2.5655077  6.327544 15.89231 6.705054e-05 2.167551e-02
# SoxN        1.4281974  7.484610 15.88021 6.748078e-05 2.167551e-02
# Tep4        1.7510950  6.435289 15.81634 6.979715e-05 2.167551e-02
# Tbp         2.1967579  5.687018 15.62166 7.736321e-05 2.288109e-02
# TotA        1.8441029  7.278051 15.37989 8.791894e-05 2.482112e-02
# CG10289    -0.8016341 10.521044 15.27596 9.289114e-05 2.508465e-02
# Root       -2.5900866  5.613942 15.12672 1.005302e-04 2.601639e-02
# cu          1.3224490  7.269899 15.00367 1.073023e-04 2.665818e-02
# PGRP-SA     2.1555869  6.415707 14.77509 1.211249e-04 2.893488e-02
# Sp7         1.9447756  5.886763 14.59283 1.334209e-04 3.069175e-02
# CG32365     2.9102548  6.957829 14.04454 1.785314e-04 3.960209e-02
# SPE         1.5362648  6.425691 13.94105 1.886340e-04 4.040020e-02
# Mtk         1.7992252  8.606391 13.60678 2.253702e-04 4.665914e-02
# Pgi         1.0033033  8.779037 13.48641 2.402981e-04 4.814489e-02
# cer        -1.1121659  7.862132 13.30139 2.652091e-04 5.147542e-02
summary(decideTests(lrt)) # 6 down 25 up
rownames(topTags(lrt, n = 35,p.value = 0.05)) -> edgeR.DE.gene
# 5 down-regulated genes & 23 up-regulated genes coming out of 
#the likelihood ratio test edgeR analysis from Starved - Fed comparison
# write.csv(as.data.frame(topTags(lrt,n = 600)),
#           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk4_8Starved-15Fed_edgeR_20201013.csv")

data.table(gene = rownames(as.data.frame(topTags(lrt,n = 28))),
           as.data.frame(topTags(lrt,n = 28)),
           method = 'edgeR.lrt', dataset = 'Lk.4' ) -> lk4.edgeR.de.tab
all.edgeR.de.tab
all.DESeq2.de.tab

# DESeq2
g.colData.df
str(g.colData.df)
as.factor(g.colData.df[,"Train_protocol"]) -> g.colData.df$Train_protocol
dds <- DESeqDataSetFromMatrix(countData = g.count.df, 
                              colData = g.colData.df, 
                              design = ~ Train_protocol)
# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
filter <- rowSums(assay(dds)>1)>5
table(filter) # FALSE 10984, TRUE 6127
# Identify the 2000 most variable genes
assay(dds.preF) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE)
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
dds.preF.analysis <- DESeq(dds.preF)
res5.S2F <- results(dds.preF.analysis, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
summary(res5.S2F)
# 12 up-regulated and 7 down-regulated DE genes (padj < 0.05)
res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
res5.S2F.Ordered
write.csv(as.data.frame(res5.S2F.Ordered),
          file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/2021.rev/lk4_8Starved-15Fed_DESeq2_20210504.csv")
res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])
S2F.gene

data.table(gene = rownames(as.data.frame(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])),
           as.data.frame(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,]),
           method = 'DESeq2', dataset = 'Lk.4' ) -> lk4.DESeq2.de.tab
all.edgeR.de.tab
all.DESeq2.de.tab
lk4.edgeR.de.tab

####
# Make a bar chart for the DE genes
g.count.df
g.colData.df
bar.order = c(rownames(g.colData.df[g.colData.df$Train_protocol=='Fed',]), rownames(g.colData.df[g.colData.df$Train_protocol=='Starved',]) )
edgeR.DE.gene
unique(c(edgeR.DE.gene,S2F.gene)) -> DE.gene #34 DE genes
lk_tpm.df[rownames(lk_tpm.df)%in%DE.gene , rownames(g.colData.df)] -> gene.matrix
lk_tpm.df['CG33926',]
str(gene.matrix)
length(DE.gene)
nrow(gene.matrix)
for (i in 1:nrow(gene.matrix)){
  gather(gene.matrix[DE.gene[i],]) -> df
  pdf(paste("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk4.only/barplot.",DE.gene[i],".23lhlk.20201012.pdf", sep = "") )
  print(ggplot(df, aes(x = key, y=value)) + geom_bar(stat = "identity", width = 0.75) + 
          theme_gray(base_size = 10) + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1 )) +
          labs(title = "15 Fed vs 8 Starved",x ="sample", y = paste("tpm of ",DE.gene[i], sep = "")) + geom_vline(xintercept=c(15.5), linetype="dashed", color = "red") +
          scale_x_discrete(limits=bar.order) #+
        #geom_text(x = 3, y = mx ,label = "Fed", color="light brown", size = 9, family="Courier", fontface="plain")
        # geom_text(x = 11, y = max(df)*3/4 ,label = "Starved", color="dark green", size = 9, family="Courier", fontface="plain")
  )
  dev.off()
}

# Make a big table including all the DE analysis results
# df1[df2, alpha := i.alpha, on = c(lsr="li", ppr="pro")]
# setkey(df1, lsr, ppr)
# setkey(df2, li, pro)
# df1[df2, alpha := i.alpha]
# setkey(all.edgeR.de.tab, gene, dataset, method)
# setkey(all.DESeq2.de.tab, gene, dataset, method)
# setkey(lk4.edgeR.de.tab, gene, dataset, method)
# setkey(lk4.DESeq2.de.tab, gene, dataset, method)
# merge(lk4.edgeR.de.tab,lk4.DESeq2.de.tab, by = c('gene'),all = T)
merge(lk4.edgeR.de.tab,lk4.DESeq2.de.tab, by = c('gene'), all = T, suffixes = c("_4.e", "_4.d") ) -> lk4.de.tab
lk4.de.tab[is.na(method_4.e),method_4.e:=method_4.d]
lk4.de.tab[method_4.d!=method_4.e,method:="edgeR.lrt,DESeq2"]
lk4.de.tab[is.na(method), method:=method_4.e]
lk4.de.tab[is.na(dataset_4.d),dataset_4.d:=dataset_4.e]
names(lk4.de.tab)[16] <- c("dataset")
lk4.de.tab[,c("method_4.e","dataset_4.e","method_4.d"):=NULL]

merge(all.edgeR.de.tab,all.DESeq2.de.tab, by = c('gene'), all = T, suffixes = c("_a.e", "_a.d") ) -> all.de.tab
all.de.tab[is.na(method_a.e),method_a.e:=method_a.d]
all.de.tab[method_a.d!=method_a.e,method:="edgeR.lrt,DESeq2"]
all.de.tab[is.na(method), method:=method_a.e]
all.de.tab[is.na(dataset_a.d),dataset_a.d:=dataset_a.e]
names(all.de.tab)[16] <- c("dataset")
all.de.tab[,c("method_a.e","dataset_a.e","method_a.d"):=NULL]

merge(lk4.de.tab, all.de.tab, by = c('gene'), all = T, suffixes = c("_lk4","_all")) -> lhlk.de.tab
# fwrite(lhlk.de.tab,"/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/20201015.lhlk.de.tab.csv",append = F)
lhlk.de.tab[is.na(dataset_all), dataset_all:=dataset_lk4]
lhlk.de.tab[dataset_all!=dataset_lk4, from.dataset:="Lk4.alone&all"]
lhlk.de.tab[is.na(from.dataset), from.dataset:=dataset_all]
lhlk.de.tab[from.dataset=="Lk.4", from.dataset:="Lk4.alone"]
lhlk.de.tab[from.dataset=="Lk.1,Lk.3,Lk.4", from.dataset:="all"]
lhlk.de.tab[,c(grep('baseMean',names(lhlk.de.tab) ),grep('logCPM',names(lhlk.de.tab) ),grep('LR',names(lhlk.de.tab) ),
            grep('PValue',names(lhlk.de.tab) ,ignore.case = T),grep('lfcSE',names(lhlk.de.tab) ),grep('stat',names(lhlk.de.tab) ),
            grep('dataset_',names(lhlk.de.tab) ) ):=NULL]

setcolorder(lhlk.de.tab,c("gene","from.dataset","method_lk4","method_all","logFC_lk4","FDR_lk4","log2FoldChange_lk4","padj_lk4",
                          "logFC_all","FDR_all","log2FoldChange_all","padj_all"))
c("gene","from.dataset","method_lk4","method_all","edgeR.log2FC_lk4","edgeR.FDR_lk4","DESeq2.log2FC_lk4","DESeq2.padj_lk4",
  "edgeR.log2FC_all","edgeR.FDR_all","DESeq2.log2FC_all","DESeq2.padj_all") -> names(lhlk.de.tab)
lhlk.de.tab
# fwrite(lhlk.de.tab,"/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/20201015.lhlk.de.tab_curated.csv",append = F)

# setdiff(all.gene_24f20s.de.dt[lhlk.de.tab, on='gene']$gene, lhlk.de.tab$gene)
all.gene_24f20s.de.dt[,c("baseMean","lfcSE","stat","pvalue","logCPM","LR","PValue"):=NULL]
all.gene_24f20s.de.dt[lhlk.de.tab, on='gene'] -> full.lhlk.de.tab
full.lhlk.de.tab[,edgeR.log2FC_all:=logFC]
full.lhlk.de.tab[,edgeR.FDR_all:=FDR]
full.lhlk.de.tab[,DESeq2.log2FC_all:=log2FoldChange]
full.lhlk.de.tab[,DESeq2.padj_all:=padj]
full.lhlk.de.tab[,c(names(full.lhlk.de.tab)[2:9]):=NULL]
full.lhlk.de.tab[is.na(DESeq2.padj_all) & !is.na(DESeq2.log2FC_all),gene] -> temp.gene
full.lhlk.de.tab[is.na(DESeq2.padj_all) & !is.na(DESeq2.log2FC_all),DESeq2.padj_all:=1]
# full.lhlk.de.tab[DESeq2.padj_all==1, gene] == temp.gene
full.lhlk.de.tab
# fwrite(full.lhlk.de.tab,"/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/20201118.lhlk.de.tab_curated.csv",append = F)

# # Replace NA by ""
# f_dowle2mx = function(DT) {
#   for (i in names(DT))
#     DT[is.na(get(i)), (i):=""]
# }
# f_dowle2mx(lk4.de.tab)
lk4.de.tab
# dt[ , new := do.call(paste, c(.SD, sep = ":"))]
# dt[, new:=do.call(paste0,.SD), .SDcols=-1]
# lk4.de.tab[, method:=do.call(paste,c(.SD, sep = ",") ), .SDcols=c("method_4.e","method_4.d")]
merge(all.edgeR.de.tab,all.DESeq2.de.tab, by = c('gene'), all = T, suffixes = c("_a.e", "_a.d") ) -> all.de.tab

all.edgeR.de.tab
all.DESeq2.de.tab
lk4.edgeR.de.tab









### Split the 44 gLHLKs into early and late groups ## ANCHOR
### 20201020  ###
### Quality control # ANCHOR
## Lk.1 + Lk.3 + Lk.4 with filter of 
#colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Assigned.reads >250000 & colData.df$Cell_id!=45

data.frame(row.names = rownames(g.colData.df), mutate(g.colData.df, time.lag = Time_anesthetized.h + Time_on_stage.h) )-> g.colData.df
median(g.colData.df$time.lag) # 3.275
g.colData.df[g.colData.df$time.lag < 3.275, "Train_protocol"] %>% table() # 11 Fed & 11 Starved
g.colData.df[g.colData.df$time.lag > 3.275, "Train_protocol"] %>% table() # 13 Fed & 9 Starved
data.frame(row.names = rownames(g.colData.df), mutate(g.colData.df, time.lag.group = "E") ) -> g.colData.df
g.colData.df[g.colData.df$time.lag > 3.275, "time.lag.group"]<- c("L")
g.colData.df[g.colData.df$time.lag.group=="E",] -> E.g.colData.df
g.colData.df[g.colData.df$time.lag.group=="L",] -> L.g.colData.df
ggplot(E.g.colData.df, aes(x = Train_protocol, y = time.lag ) )+ geom_point(size = 2) + theme_gray(base_size = 14) + 
  labs(title="Early subset (<3.275 hr)", x ="State", y = "Time lag (hr)") 
ggplot(L.g.colData.df, aes(x = Train_protocol, y = time.lag ) )+ geom_point(size = 2) + theme_gray(base_size = 14) + 
  labs(title="Delayed subset (>3.275 hr)", x ="State", y = "Time lag (hr)") 

# Early subset
## DE analysis
# edgeR
## Test whether CEL-Seq data is zero-inflated using Scatterplots of 
##   the estimated biological coefficient of variation (BCV, defined 
##   as the square root of the negative binomial dispersion parameter φ)
##   against average log counts per million (CPM) computed using edgeR
nrow(E.g.colData.df)
g.count.df[,rownames(E.g.colData.df)] -> E.g.count.df
ncol(E.g.count.df)

c(as.vector(E.g.colData.df$Train_protocol) ) -> group
y <- DGEList(counts=E.g.count.df, group=group)
y$samples
y$counts['coGFP',]

# Filtering
keep <- rowSums(cpm(y)>1) >= 11 # Total 22 samples, use 11 as cutoff
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)
# 6310 genes passed

# Normalization for RNA composition
y <- calcNormFactors(y)
y$samples

# group lib.size norm.factors
# lk.1_043 Starved    38512    1.0949150
# lk.1_044 Starved    64193    0.9772050
# lk.1_046     Fed    65187    0.9452060
# lk.1_053 Starved    68287    0.8683499
# lk.1_054     Fed    45789    0.8919409
# lk.1_055     Fed    24416    1.0776399
# lk.1_060 Starved    21555    1.0729303
# lk.3_085 Starved    31811    1.0882191
# lk.3_086 Starved    47122    0.9707269
# lk.3_096     Fed    40432    1.0553486
# lk.4_098     Fed    32724    1.0046352
# lk.4_099     Fed   118383    1.0222073
# lk.4_103 Starved   118844    1.0358514
# lk.4_104 Starved    46370    0.9151429
# lk.4_105 Starved    37482    1.0280011
# lk.4_107 Starved   145413    1.0314801
# lk.4_108 Starved    75413    0.9465041
# lk.4_117     Fed    37098    1.0299698
# lk.4_118     Fed    38282    1.0701447
# lk.4_119     Fed    42073    1.0077284
# lk.4_121     Fed    47935    0.9134807
# lk.4_122     Fed    51969    0.9990827

# plotBCV after estimating common dispersion and tagwise dispersions
#  in one run (recommended):
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/BCV.plot_22.E.gLHLK_Train.Prot_20201020.pdf")
# plotBCV(y)
# dev.off()

et <- exactTest(y, pair = c("Fed","Starved"))
topTags(et, n = 20)
summary(decideTests(et)) # at 5% FDR
# 3 down & 2 up genes coming out of default edgeR analysis from Starved - Fed comparison
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
topTags(lrt, n = 35)
summary(decideTests(lrt))
# Coefficient:  groupStarved 
# logFC    logCPM        LR       PValue          FDR
# fit       -2.240486  7.891263 32.527068 1.175439e-08 7.417022e-05
# hll        3.677558  6.483294 30.322645 3.658300e-08 1.154194e-04
# Got2       1.837062  8.996505 22.586902 2.008493e-06 4.224530e-03
# CG5773    -2.117109  7.316507 21.118257 4.317946e-06 6.811560e-03
# Yp2       -1.146775 12.363764 20.646880 5.522710e-06 6.969660e-03
# CG31075    1.860315  6.850599 14.093352 1.739575e-04 1.829453e-01
rownames(topTags(lrt, n = 35,p.value = 0.05)) -> E.edgeR.DE.gene
# 3 down-regulated genes & 2 up-regulated genes coming out of 
#the likelihood ratio test edgeR analysis from Starved - Fed comparison
write.csv(as.data.frame(topTags(lrt,n = 250)),
          file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk_early.11Starved-11Fed_edgeR_20201020.csv")
data.table(gene = rownames(as.data.frame(topTags(lrt,n = 5))),as.data.frame(topTags(lrt,n = 5)),
           method = 'edgeR.lrt', dataset = 'Lk.1,Lk.3,Lk.4', stage = 'early' ) -> E.all.edgeR.de.tab

# DESeq2
E.g.colData.df
str(E.g.colData.df)
as.factor(E.g.colData.df[,"Train_protocol"]) -> E.g.colData.df$Train_protocol
dds <- DESeqDataSetFromMatrix(countData = E.g.count.df, 
                              colData = E.g.colData.df, 
                              design = ~ Train_protocol)
# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
filter <- rowSums(assay(dds)>1)>5
table(filter) # FALSE 7486, TRUE 5861
# Identify the 2000 most variable genes
assay(dds.preF) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE)
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
dds.preF.analysis <- DESeq(dds.preF)
res5.S2F <- results(dds.preF.analysis, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
summary(res5.S2F)
# 1 up-regulated and 3 down-regulated DE genes (padj < 0.05)
res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
res5.S2F.Ordered
write.csv(as.data.frame(res5.S2F.Ordered),
          file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk_E.11Starved-11Fed_DESeq2_20201020.csv")
res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])
S2F.gene
data.table(gene = rownames(as.data.frame(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])),
           as.data.frame(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,]),
           method = 'DESeq2', dataset = 'Lk.1,Lk.3,Lk.4', stage = 'early' ) -> E.all.DESeq2.de.tab
E.all.edgeR.de.tab
union(E.all.DESeq2.de.tab$gene, E.all.edgeR.de.tab$gene)
# Delayed stage
## DE analysis
# edgeR
## Test whether CEL-Seq data is zero-inflated using Scatterplots of 
##   the estimated biological coefficient of variation (BCV, defined 
##   as the square root of the negative binomial dispersion parameter φ)
##   against average log counts per million (CPM) computed using edgeR
nrow(L.g.colData.df)
g.count.df[,rownames(L.g.colData.df)] -> L.g.count.df
ncol(L.g.count.df)

c(as.vector(L.g.colData.df$Train_protocol) ) -> group
y <- DGEList(counts=L.g.count.df, group=group)
y$samples
y$counts['coGFP',]

# Filtering
keep <- rowSums(cpm(y)>1) >= 11 # Total 22 samples, use 11 as cutoff
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)
# 6422 genes passed

# Normalization for RNA composition
y <- calcNormFactors(y)
y$samples

# group lib.size norm.factors
# lk.1_049 Starved    33214    1.0603736
# lk.1_050 Starved    51754    0.8905533
# lk.1_052 Starved    60716    1.0249134
# lk.1_057     Fed    59411    1.0243453
# lk.1_058     Fed    77840    1.0167920
# lk.3_078     Fed    33928    1.0553097
# lk.3_088     Fed    51572    0.9442309
# lk.3_089     Fed    36085    1.0724915
# lk.3_092 Starved    60876    0.8756820
# lk.3_100 Starved    38419    1.1136844
# lk.3_101 Starved    48308    0.8915592
# lk.4_028 Starved    85161    1.0228459
# lk.4_029     Fed    32065    1.1615188
# lk.4_102 Starved    42698    0.8649191
# lk.4_110     Fed    65621    1.0127130
# lk.4_111     Fed   101016    1.0638842
# lk.4_112     Fed    79146    0.9788679
# lk.4_113     Fed    56187    1.0508007
# lk.4_114     Fed    71539    0.9699215
# lk.4_115     Fed    36623    1.0156219
# lk.4_116     Fed    36986    0.9777725
# lk.4_123 Starved    94459    0.9733324

# plotBCV after estimating common dispersion and tagwise dispersions
#  in one run (recommended):
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)
pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/BCV.plot_22.L.gLHLK_Train.Prot_20201020.pdf")
plotBCV(y)
dev.off()

et <- exactTest(y, pair = c("Fed","Starved"))
topTags(et, n = 20)
summary(decideTests(et)) # at 5% FDR
# 0 down & 0 up genes coming out of default edgeR analysis from Starved - Fed comparison
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
topTags(lrt, n = 35)
summary(decideTests(lrt))

rownames(topTags(lrt, n = 35,p.value = 0.05)) -> L.edgeR.DE.gene
# Nothing coming out of 
#the likelihood ratio test edgeR analysis from Starved - Fed comparison
# write.csv(as.data.frame(topTags(lrt,n = 100)),
#           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk_early.11Starved-11Fed_edgeR_20201020.csv")
data.table(gene = rownames(as.data.frame(topTags(lrt,n = 0))),as.data.frame(topTags(lrt,n = 0)),
           method = 'edgeR.lrt', dataset = 'Lk.1,Lk.3,Lk.4', stage = 'delayed' ) -> L.all.edgeR.de.tab

















### plus the stringent counted read num fitler (40k) # ANCHOR
# NO DE gene from edgeR nor DESeq2
(median(colData.df$Counted.reads) + min(colData.df$Counted.reads) )/2 # 22843.5
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Counted.reads >40000, ] -> g.colData.df # In DPM project, 0.035 and 0.8 as thresholds
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Counted.reads >40000
           , "Train_protocol"] %>% table()
# 6 Fed vs 10 Starved 
rownames(setdiff(colData.df,g.colData.df)) 
# [1] "lk.1_041" "lk.1_042" "lk.1_045" "lk.1_047" "lk.1_049" "lk.1_051" "lk.1_055" "lk.1_056" "lk.1_060" "lk.1_061" "lk.3_067"
# [12] "lk.3_078" "lk.3_084" "lk.3_085" "lk.3_089" "lk.3_090" "lk.3_091" "lk.3_093" "lk.3_094" "lk.3_095" "lk.3_099"

colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Counted.reads >40000, 
           "detected.gene.count"] %>% mean() # 6389.5
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Counted.reads >40000, 
           "detected.gene.count"] %>% median() # 6214.5
rownames(g.colData.df) -> g.lk.id

data.frame(lk_counts[,2:38], row.names = lk_counts$gene_id) -> count.df
count.df[,g.lk.id] -> g.count.df
ncol(g.count.df)


pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/plot_bar_16.lk_de.g.20200922.pdf", width = 9)
ggplot(data = g.colData.df, aes(x = rownames(g.colData.df), y=detected.gene.count) ) + geom_bar(stat = "identity", width = 0.6) +
  scale_x_discrete(limits=rownames(g.colData.df) ) + theme_gray(base_size = 18) +
  theme(axis.text.x = element_text(angle=45 , hjust = 1, vjust = 1) ) +
  labs(title="6 Fed vs 10 Starved LHLKs", x ="Cell_id", y = "Detected gene count") + geom_hline(yintercept=median(as.vector(colSums(g.count.df>0))), linetype="dashed", color = "red") +
  geom_text(x = 13, y = 7500,label = paste('median = ',median(as.vector(colSums(g.count.df>0))), sep = ''), color="red", size = 8, family="Courier", fontface="plain")
dev.off()

# edgeR
## Test whether CEL-Seq data is zero-inflated using Scatterplots of 
##   the estimated biological coefficient of variation (BCV, defined 
##   as the square root of the negative binomial dispersion parameter φ)
##   against average log counts per million (CPM) computed using edgeR
nrow(g.colData.df)
ncol(g.count.df)

c(as.vector(g.colData.df$Train_protocol) ) -> group
y <- DGEList(counts=g.count.df, group=group)
y$samples
y$counts['coGFP',]

# Filtering
keep <- rowSums(cpm(y)>1) >= 8 # Total 16 samples, use 8 as cutoff
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)
# 6487 genes passed

# Normalization for RNA composition
y <- calcNormFactors(y)
y$samples

# group lib.size norm.factors
# lk.1_043 Starved    38593    1.0089035
# lk.1_044 Starved    64390    1.0932431
# lk.1_046     Fed    65439    1.0369481
# lk.1_050 Starved    51856    0.9505006
# lk.1_052 Starved    60776    1.1191393
# lk.1_053 Starved    68573    0.9074231
# lk.1_054     Fed    45971    0.9011570
# lk.1_057     Fed    59551    1.1131542
# lk.1_058     Fed    77953    0.9716476
# lk.3_086 Starved    47238    1.0565165
# lk.3_087 Starved    41684    0.9311858
# lk.3_088     Fed    51603    0.9500198
# lk.3_092 Starved    60919    0.9460795
# lk.3_096     Fed    40565    1.0418508
# lk.3_100 Starved    38425    1.0330183
# lk.3_101 Starved    48266    0.9765156

# plotBCV after estimating common dispersion and tagwise dispersions
#  in one run (recommended):
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)
pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/BCV.plot_16gLHLK_Train.Prot_20200922.pdf")
# y <- estimateCommonDisp(y)
# y <- estimateTrendedDisp(y)
# y <- estimateTagwiseDisp(y)
plotBCV(y)
dev.off()

et <- exactTest(y, pair = c("Fed","Starved"))
topTags(et, n = 20)
summary(decideTests(et)) # at 5% FDR

# two down genes coming out of default edgeR analysis from Starved - Fed comparison
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
topTags(lrt, n = 30)
summary(decideTests(lrt))
# Coefficient:  groupStarved 
# logFC   logCPM        LR       PValue       FDR
# CG5773      -2.048353 7.512406 18.334640 1.853069e-05 0.1185494
# fit         -2.056829 7.992012 16.734571 4.299045e-05 0.1185494
# CG5500      -2.428844 6.058355 16.273521 5.482477e-05 0.1185494
# CG3699       3.965013 6.835021 15.683878 7.485956e-05 0.1214035
# CG1092      -2.603787 5.874492 14.012958 1.815551e-04 0.2355496
rownames(topTags(lrt, n = 30,p.value = 0.05)) -> edgeR.DE.gene
# No gene coming out of the likelihood ratio test edgeR analysis from Starved - Fed comparison
# write.csv(as.data.frame(topTags(lrt, n = 30)),
#           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk_10Starved-6Fed_edgeR_20200922.csv")


# DESeq2
g.colData.df
str(g.colData.df)
as.factor(g.colData.df[,"Train_protocol"]) -> g.colData.df$Train_protocol
dds <- DESeqDataSetFromMatrix(countData = g.count.df, 
                              colData = g.colData.df, 
                              design = ~ Train_protocol)
# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
filter <- rowSums(assay(dds)>1)>5
table(filter) # FALSE 6295, TRUE 5093
# Identify the 2000 most variable genes
assay(dds.preF) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE)
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
dds.preF.analysis <- DESeq(dds.preF)
res5.S2F <- results(dds.preF.analysis, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
summary(res5.S2F)
res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
res5.S2F.Ordered
# write.csv(as.data.frame(res5.S2F.Ordered),
#           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk_10Starved-6Fed_DESeq2_20200922.csv")
res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])
S2F.gene


### Loosen filter plus the counted read num inside lk.3 # ANCHOR
# no gene coming out of either pipelines # 20200921
(median(colData.df$Counted.reads) + min(colData.df$Counted.reads) )/2 # 22843.5
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Counted.reads >20000 & 
             colData.df$Library=='Lk.3', ] -> g.colData.df 
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Counted.reads >20000 & 
             colData.df$Library=='Lk.3', "Train_protocol"] %>% table()
# 6 Fed vs 9 Starved 
rownames(setdiff(colData.df,g.colData.df)) # Remove  "lk.3_067" "lk.3_090" "lk.3_093"

colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Counted.reads >20000 & 
             colData.df$Library=='Lk.3', "detected.gene.count"] %>% mean() # 5550.733
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.80 & colData.df$Counted.reads >20000 & 
             colData.df$Library=='Lk.3', "detected.gene.count"] %>% median() # 5612
rownames(g.colData.df) -> g.lk.id
data.frame(lk_counts[,2:38], row.names = lk_counts$gene_id) -> count.df
count.df[,g.lk.id] -> g.count.df
ncol(g.count.df)


pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/plot_bar_6F-9S.lk3_de.g.20200921.pdf", width = 10)
ggplot(data = g.colData.df, aes(x = rownames(g.colData.df), y=detected.gene.count) ) + geom_bar(stat = "identity", width = 0.6) +
  scale_x_discrete(limits=rownames(g.colData.df) ) + theme_gray(base_size = 18) +
  theme(axis.text.x = element_text(angle=45 , hjust = 1, vjust = 1) ) +
  labs(title="6 Fed vs 9 Starved LHLKs", x ="Cell_id", y = "Detected gene count") + geom_hline(yintercept=median(as.vector(colSums(g.count.df>0))), linetype="dashed", color = "red") +
  geom_text(x = 19, y = 7500,label = paste('median = ',median(as.vector(colSums(g.count.df>0))), sep = ''), color="red", size = 8, family="Courier", fontface="plain")
dev.off()

# edgeR
## Test whether CEL-Seq data is zero-inflated using Scatterplots of 
##   the estimated biological coefficient of variation (BCV, defined 
##   as the square root of the negative binomial dispersion parameter φ)
##   against average log counts per million (CPM) computed using edgeR
nrow(g.colData.df)
ncol(g.count.df)

c(as.vector(g.colData.df$Train_protocol) ) -> group
y <- DGEList(counts=g.count.df, group=group)
y$samples
y$counts['coGFP',]

# Filtering
keep <- rowSums(cpm(y)>1) >= 8 # Total 15 samples, use 8 as cutoff
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)
# 5260 genes passed

# Normalization for RNA composition
y <- calcNormFactors(y)
y$samples

# group lib.size norm.factors
# lk.3_078     Fed    33200    0.9304337
# lk.3_084 Starved    25222    1.0257455
# lk.3_085 Starved    31277    0.9857463
# lk.3_086 Starved    46281    1.0783683
# lk.3_087 Starved    40803    0.8958018
# lk.3_088     Fed    50515    0.8920454
# lk.3_089     Fed    35516    1.0233552
# lk.3_091 Starved    27273    0.9803504
# lk.3_092 Starved    59660    0.8971209
# lk.3_094     Fed    18779    1.0936449
# lk.3_095     Fed    23770    1.1234161
# lk.3_096     Fed    39606    1.0634774
# lk.3_099 Starved    21338    1.0444142
# lk.3_100 Starved    37660    1.0424492
# lk.3_101 Starved    47453    0.9634122

# plotBCV after estimating common dispersion and tagwise dispersions
#  in one run (recommended):
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)
pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/BCV.plot_15gLHLK.Lk3_Train.Prot_20200921.pdf")
# y <- estimateCommonDisp(y)
# y <- estimateTrendedDisp(y)
# y <- estimateTagwiseDisp(y)
plotBCV(y)
dev.off()

et <- exactTest(y, pair = c("Fed","Starved"))
topTags(et, n = 20)
summary(decideTests(et)) # at 5% FDR
# No gene coming out of default edgeR analysis from Starved - Fed comparison

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
topTags(lrt, n = 30)
summary(decideTests(lrt))
# No gene coming out of the likelihood ratio test edgeR analysis from Starved - Fed comparison

# DESeq2
g.colData.df
str(g.colData.df)
as.factor(g.colData.df[,"Train_protocol"]) -> g.colData.df$Train_protocol
dds <- DESeqDataSetFromMatrix(countData = g.count.df, 
                              colData = g.colData.df, 
                              design = ~ Train_protocol)
# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
filter <- rowSums(assay(dds)>1)>5
table(filter) # FALSE 7433, TRUE 3955
# Identify the 2000 most variable genes
assay(dds.preF) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE)
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
dds.preF.analysis <- DESeq(dds.preF)
res5.S2F <- results(dds.preF.analysis, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
summary(res5.S2F)
res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
res5.S2F.Ordered
# No gene coming out of the pipeline


#### Succeeded cutoff (0.09, 0.85) but I have changed my mind after getting to know that CEL-Seq paper didn't take # ANCHOR
# the ERCCs smaller than 10 counts for the linearity calculation. The increased technical noise is to blame.
# 3 down (CG5773, CG1092, fit) & 3 up (hll, Got2, CG3699) from edgeR
# 3 down from DESeq2 (CG5773, Hml, fit)

colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.85, ] -> g.colData.df # In DPM project, 0.035 and 0.8 as thresholds
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.85, "Train_protocol"] %>% table()
# 12 Fed vs 14 Starved # Remove 041 (MBON), 042, 047, 051, 060 and 061
rownames(setdiff(colData.df,g.colData.df)) # Remove 041 (MBON),"lk.1_042" "lk.1_047" "lk.1_051" 
#"lk.1_060" "lk.1_061" "lk.3_067" "lk.3_086" "lk.3_093" "lk.3_095" "lk.3_096"

colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.85, "detected.gene.count"] %>% mean() # 5707.077
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.85, "detected.gene.count"] %>% median() # 5665
rownames(g.colData.df) -> g.lk.id
data.frame(lk_counts[,2:38], row.names = lk_counts$gene_id) -> count.df
count.df[,g.lk.id] -> g.count.df
ncol(g.count.df)


pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/plot_bar_26.lk_de.g.20200921.pdf", width = 10)
ggplot(data = g.colData.df, aes(x = rownames(g.colData.df), y=detected.gene.count) ) + geom_bar(stat = "identity", width = 0.6) +
  scale_x_discrete(limits=rownames(g.colData.df) ) + theme_gray(base_size = 18) +
  theme(axis.text.x = element_text(angle=45 , hjust = 1, vjust = 1) ) +
  labs(title="12 Fed vs 14 Starved LHLKs", x ="Cell_id", y = "Detected gene count") + geom_hline(yintercept=median(as.vector(colSums(g.count.df>0))), linetype="dashed", color = "red") +
  geom_text(x = 19, y = 7500,label = paste('median = ',median(as.vector(colSums(g.count.df>0))), sep = ''), color="red", size = 8, family="Courier", fontface="plain")
dev.off()

# edgeR
## Test whether CEL-Seq data is zero-inflated using Scatterplots of 
##   the estimated biological coefficient of variation (BCV, defined 
##   as the square root of the negative binomial dispersion parameter φ)
##   against average log counts per million (CPM) computed using edgeR
nrow(g.colData.df)
ncol(g.count.df)

c(as.vector(g.colData.df$Train_protocol) ) -> group
y <- DGEList(counts=g.count.df, group=group)
y$samples
y$counts['coGFP',]

# Filtering
keep <- rowSums(cpm(y)>1) >= 13 # Total 26 samples, use 13 as cutoff
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)
# 5612 genes passed

# Normalization for RNA composition
y <- calcNormFactors(y)
y$samples

# group lib.size norm.factors
# lk.1_043 Starved    37932    1.1100937
# lk.1_044 Starved    63528    1.0051432
# lk.1_045     Fed    24558    0.9566190
# lk.1_046     Fed    64403    0.9689408
# lk.1_049 Starved    32786    1.0488998
# lk.1_050 Starved    51113    0.8877682
# lk.1_052 Starved    59881    1.0683671
# lk.1_053 Starved    67659    0.8916649
# lk.1_054     Fed    45330    0.9098947
# lk.1_055     Fed    24095    1.0751306
# lk.1_056     Fed    23690    1.1444603
# lk.1_057     Fed    58484    1.0599479
# lk.1_058     Fed    76771    0.9866806
# lk.3_078     Fed    33451    1.0191302
# lk.3_084 Starved    25327    1.0159486
# lk.3_085 Starved    31426    1.0861189
# lk.3_087 Starved    41074    0.9183901
# lk.3_088     Fed    50862    0.8704779
# lk.3_089     Fed    35694    1.0333668
# lk.3_090     Fed    17336    1.1050230
# lk.3_091 Starved    27188    1.0260363
# lk.3_092 Starved    60112    0.9151461
# lk.3_094     Fed    18829    1.0229121
# lk.3_099 Starved    21394    0.9773413
# lk.3_100 Starved    37888    1.0567811
# lk.3_101 Starved    47644    0.9141004

# plotBCV after estimating common dispersion and tagwise dispersions
#  in one run (recommended):
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)
pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/BCV.plot_26gLHLK_Train.Prot_20200921.pdf")
# y <- estimateCommonDisp(y)
# y <- estimateTrendedDisp(y)
# y <- estimateTagwiseDisp(y)
plotBCV(y)
dev.off()

et <- exactTest(y, pair = c("Fed","Starved"))
topTags(et, n = 20)
summary(decideTests(et)) # at 5% FDR
# Comparison of groups:  Starved-Fed 
# logFC    logCPM       PValue          FDR
# CG5773  -1.9873517  7.828933 3.440671e-08 0.0001930905
# CG1092  -3.1404349  6.356220 2.779196e-07 0.0007798423
# fit     -1.8546323  8.340499 2.123260e-05 0.0397191151
# Got2     1.5597227  9.125555 3.755731e-05 0.0526929067
# hll      2.5936583  6.721367 5.290175e-05 0.0593769228
# CG3699   2.7329454  6.945088 9.088265e-05 0.0850055763
# Yp2     -1.0923308 12.578354 1.385461e-04 0.1110743857
# Three genes coming out of default edgeR analysis from Starved - Fed comparison
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
topTags(lrt, n = 30)
summary(decideTests(lrt))
# Coefficient:  groupStarved 
# logFC    logCPM        LR       PValue          FDR
# CG5773      -1.9873517  7.828933 31.191027 2.338445e-08 0.0001312335
# CG1092      -3.1404349  6.356220 27.682017 1.429858e-07 0.0004012181
# fit         -1.8546323  8.340499 18.681026 1.545124e-05 0.0289041248
# hll          2.5936583  6.721367 17.506792 2.862832e-05 0.0329330234
# Got2         1.5597227  9.125555 17.301448 3.189438e-05 0.0329330234
# CG3699       2.7329454  6.945088 17.113565 3.520993e-05 0.0329330234
# Yp2         -1.0923308 12.578354 14.563945 1.354822e-04 0.1086180340
rownames(topTags(lrt, n = 30,p.value = 0.05)) -> edgeR.DE.gene
# Three down-regulated genes (CG5773, CG1092, fit) & three up-regulated genes (hll, Got2, CG3699) coming out of 
#the likelihood ratio test edgeR analysis from Starved - Fed comparison


# DESeq2
g.colData.df
str(g.colData.df)
as.factor(g.colData.df[,"Train_protocol"]) -> g.colData.df$Train_protocol
dds <- DESeqDataSetFromMatrix(countData = g.count.df, 
                              colData = g.colData.df, 
                              design = ~ Train_protocol)
# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
filter <- rowSums(assay(dds)>1)>5
table(filter) # FALSE 6002, TRUE 5386
# Identify the 2000 most variable genes
assay(dds.preF) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE)
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
dds.preF.analysis <- DESeq(dds.preF)
res5.S2F <- results(dds.preF.analysis, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
summary(res5.S2F)
res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
res5.S2F.Ordered
# write.csv(as.data.frame(res5.S2F.Ordered),
#           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk_14Starved-12Fed__20200921.csv")
res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])
S2F.gene




----
# put all 36 LHLKs
# DESeq2
dds <- DESeqDataSetFromMatrix(countData = count.df[,-1], 
                                colData = colData.df[-1,], 
                                design = ~ Train_protocol)
colData.df[-1, "Train_protocol"] %>% table() # 18 Fed vs 18 Starved

# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
filter <- rowSums(assay(dds)>1)>5
table(filter) # FALSE 8600, TRUE 4430
# Identify the 2000 most variable genes
assay(dds.preF) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE)
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
dds.preF.analysis <- DESeq(dds.preF)
res5.S2F <- results(dds.preF.analysis, alpha = 0.05, contrast = c("Train_protocol","Starved","Fed"))
summary(res5.S2F)
res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
res5.S2F.Ordered
# write.csv(as.data.frame(res5.S2F.Ordered),
#           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk_18Starved-18Fed_20200921.csv")
res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.10,])
S2F.gene


# Loosen cutoff (0.09, 0.8) 15 F vs 16 S
# 3 down & 4 up from edgeR
# 2 down from  DESeq2

colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.8, ] -> g.colData.df # In DPM project, 0.035 and 0.8 as thresholds
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.8, "Train_protocol"] %>% table()
# 15 Fed vs 16 Starved 
rownames(setdiff(colData.df,g.colData.df)) # Remove 041 (MBON),"lk.1_042" "lk.1_047" "lk.1_051" "lk.1_061" "lk.3_093
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.8, "detected.gene.count"] %>% mean() # 5651.71
colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.8, "detected.gene.count"] %>% median() # 5680 or 5613 after removing the ERCC
rownames(g.colData.df) -> g.lk.id
data.frame(lk_counts[,2:38], row.names = lk_counts$gene_id) -> count.df
count.df[,g.lk.id] -> g.count.df
ncol(g.count.df)


# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/plot_bar_31.lk_de.g.20200921.pdf", width = 11)
# ggplot(data = g.colData.df, aes(x = rownames(g.colData.df), y=detected.gene.count) ) + geom_bar(stat = "identity", width = 0.6) +
#   scale_x_discrete(limits=rownames(g.colData.df) ) + theme_gray(base_size = 18) +
#   theme(axis.text.x = element_text(angle=45 , hjust = 1, vjust = 1) ) +
#   labs(title="15 Fed vs 16 Starved LHLKs", x ="Cell_id", y = "Detected gene count") + geom_hline(yintercept=median(as.vector(colSums(g.count.df>0))), linetype="dashed", color = "red") +
#   geom_text(x = 24, y = 7500,label = paste('median = ',median(as.vector(colSums(g.count.df>0))), sep = ''), color="red", size = 8, family="Courier", fontface="plain")
# dev.off()


# edgeR
## Test whether CEL-Seq data is zero-inflated using Scatterplots of 
##   the estimated biological coefficient of variation (BCV, defined 
##   as the square root of the negative binomial dispersion parameter φ)
##   against average log counts per million (CPM) computed using edgeR
nrow(g.colData.df)
ncol(g.count.df)

c(as.vector(g.colData.df$Train_protocol) ) -> group
y <- DGEList(counts=g.count.df, group=group)
y$samples
y$counts['coGFP',]

# Filtering
keep <- rowSums(cpm(y)>1) >= 15 # Total 31 samples, use 15 as cutoff
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)
# 5641 genes passed

# Normalization for RNA composition
y <- calcNormFactors(y)
y$samples

# group lib.size norm.factors
# lk.1_043 Starved    37940    1.0961509
# lk.1_044 Starved    63519    0.9811497
# lk.1_045     Fed    24526    0.9386740
# lk.1_046     Fed    64420    0.9585054
# lk.1_049 Starved    32794    1.0442344
# lk.1_050 Starved    51092    0.8662504
# lk.1_052 Starved    59923    1.0448216
# lk.1_053 Starved    67632    0.8729282
# lk.1_054     Fed    45296    0.8879337
# lk.1_055     Fed    24074    1.0541183
# lk.1_056     Fed    23707    1.1209559
# lk.1_057     Fed    58504    1.0369438
# lk.1_058     Fed    76832    0.9830441
# lk.1_060 Starved    21211    1.1384601
# lk.3_067     Fed    15570    1.0596509
# lk.3_078     Fed    33453    1.0062796
# lk.3_084 Starved    25319    1.0060350
# lk.3_085 Starved    31429    1.0961693
# lk.3_086 Starved    46628    0.9886291
# lk.3_087 Starved    41072    0.8961318
# lk.3_088     Fed    50902    0.8572434
# lk.3_089     Fed    35721    1.0522420
# lk.3_090     Fed    17347    1.0841923
# lk.3_091 Starved    27214    0.9706363
# lk.3_092 Starved    60094    0.8926580
# lk.3_094     Fed    18834    1.0215578
# lk.3_095     Fed    23851    1.1905245
# lk.3_096     Fed    39892    1.0713957
# lk.3_099 Starved    21389    0.9556161
# lk.3_100 Starved    37918    1.0487541
# lk.3_101 Starved    47691    0.8924111

# plotBCV after estimating common dispersion and tagwise dispersions
#  in one run (recommended):
design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)
pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/BCV.plot_31gLHLK_Train.Prot_20200921.pdf")
# y <- estimateCommonDisp(y)
# y <- estimateTrendedDisp(y)
# y <- estimateTagwiseDisp(y)
plotBCV(y)
dev.off()

et <- exactTest(y, pair = c("Fed","Starved"))
topTags(et, n = 20)
summary(decideTests(et)) # at 5% FDR
# Comparison of groups:  Starved-Fed 
# logFC    logCPM       PValue        FDR
# CG5773   -1.8445072  7.723193 4.427580e-06 0.01527824
# hll       2.7085032  6.691846 7.509760e-06 0.01527824
# fit      -1.8687056  8.331419 8.394710e-06 0.01527824
# Got2      1.5562992  9.039628 1.083371e-05 0.01527824
# Yp2      -1.0988808 12.612168 1.803255e-05 0.01891967
# Tbh       2.1249389  7.006501 2.012374e-05 0.01891967
# CG3699    2.6472428  6.891132 4.377779e-05 0.03527865
# CG9119    1.6699397  8.264607 1.203107e-04 0.08483411
# 3 down & 4 up genes coming out of default edgeR analysis from Starved - Fed comparison
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
topTags(lrt, n = 30)
summary(decideTests(lrt))
# Coefficient:  groupStarved 
# logFC    logCPM        LR       PValue         FDR
# CG5773     -1.8445072  7.723193 21.973975 2.763727e-06 0.009187596
# hll         2.7085032  6.691846 21.658616 3.257435e-06 0.009187596
# fit        -1.8687056  8.331419 20.316005 6.564915e-06 0.012344228
# Got2        1.5562992  9.039628 19.465548 1.024306e-05 0.014445278
# Yp2        -1.0988808 12.612168 18.431263 1.761443e-05 0.017497550
# CG3699      2.6472428  6.891132 18.326390 1.861112e-05 0.017497550
# Tbh         2.1249389  7.006501 16.188035 5.735523e-05 0.046220125
# CG9119      1.6699397  8.264607 15.251154 9.411925e-05 0.066365839
# AdSL        1.6022236  6.896738 11.822234 5.852758e-04 0.320413599
rownames(topTags(lrt, n = 30,p.value = 0.05)) -> edgeR.DE.gene  

# Three down-regulated genes (CG5773, fit, Yp2) & four up-reg genes (hll, Got2, CG3699, Tbh) 
#coming out of the likelihood ratio test edgeR analysis from Starved - Fed comparison


# DESeq2
g.colData.df
str(g.colData.df)
as.factor(g.colData.df[,"Train_protocol"]) -> g.colData.df$Train_protocol
dds <- DESeqDataSetFromMatrix(countData = g.count.df, 
                              colData = g.colData.df, 
                              design = ~ Train_protocol)
# First, we filter out the lowly expressed genes, by removing those genes that 
#   do not have at least 1 reads in at least half of the samples.
dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
filter <- rowSums(assay(dds)>1)>5
table(filter) # FALSE 5750, TRUE 5638
# Identify the 2000 most variable genes
assay(dds.preF) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(dds.preF)
vars <- sort(vars, decreasing = TRUE)
head(vars)
dds.preF <- dds.preF[names(vars)[1:2000],]
dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
dds.preF.analysis <- DESeq(dds.preF)
res5.S2F <- results(dds.preF.analysis, alpha = 0.10, contrast = c("Train_protocol","Starved","Fed"))
summary(res5.S2F)
res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
res5.S2F.Ordered
# write.csv(as.data.frame(res5.S2F.Ordered),
#           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk_16Starved-15Fed__20200921.csv")
res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])
S2F.gene




####
# # additional cutoff of dtected.gene.count of 4500
# 
# colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.85 & colData.df$detected.gene.count > 4500, ] -> g.colData.df # In DPM project, 0.035 and 0.8 as thresholds
# colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.85 & colData.df$detected.gene.count > 4500, "Train_protocol"] %>% table()
# # 6 Fed vs 6 Starved
# colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.85 & colData.df$detected.gene.count > 4500, "detected.gene.count"] %>% mean() # 6043.308
# colData.df[colData.df$ERCCoverTranscript_ratio<0.09 & colData.df$Spike.in.linearity>0.85 & colData.df$detected.gene.count > 4500, "detected.gene.count"] %>% median() # 6126
# rownames(g.colData.df) -> g.lk.1.id
# data.frame(lk1.counts[,2:20], row.names = lk1.counts$gene_id) -> count.df
# count.df[,g.lk.1.id] -> g.count.df
# ncol(g.count.df)
# 
# 
# # pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/plot_bar_13.lk.1_de.g.20200702.pdf", width = 10)
# # ggplot(data = g.colData.df, aes(x = rownames(g.colData.df), y=detected.gene.count) ) + geom_bar(stat = "identity", width = 0.6) +
# #   scale_x_discrete(limits=rownames(g.colData.df) ) + theme_gray(base_size = 18) +
# #   theme(axis.text.x = element_text(angle=45 , hjust = 1, vjust = 1) ) +
# #   labs(title="7 Fed vs 6 Starved LHLKs", x ="Cell_id", y = "Detected gene count") + geom_hline(yintercept=median(as.vector(colSums(g.count.df>0))), linetype="dashed", color = "red") +
# #   geom_text(x = 7, y = 7500,label = paste('median = ',median(as.vector(colSums(g.count.df>0))), sep = ''), color="red", size = 9, family="Courier", fontface="plain")
# # dev.off()
# 
# 
# # edgeR
# ## Test whether CEL-Seq data is zero-inflated using Scatterplots of 
# ##   the estimated biological coefficient of variation (BCV, defined 
# ##   as the square root of the negative binomial dispersion parameter φ)
# ##   against average log counts per million (CPM) computed using edgeR
# nrow(g.colData.df)
# ncol(g.count.df)
# 
# c(as.vector(g.colData.df$Train_protocol) ) -> group
# y <- DGEList(counts=g.count.df, group=group)
# y$samples
# y$counts['coGFP',]
# 
# # Filtering
# keep <- rowSums(cpm(y)>1) >= 7 # Total 15 samples, use 7 as cutoff
# y <- y[keep, , keep.lib.sizes=FALSE]
# nrow(y$counts)
# # 5730 genes passed
# 
# # Normalization for RNA composition
# y <- calcNormFactors(y)
# y$samples
# 
# # group lib.size norm.factors
# # lk.1_043 Starved    38142    1.1189227
# # lk.1_044 Starved    63750    1.0348758
# # lk.1_045     Fed    24686    1.0191239
# # lk.1_046     Fed    64618    0.9006138
# # lk.1_049 Starved    33095    1.0370791
# # lk.1_050 Starved    51175    0.9889356
# # lk.1_052 Starved    60161    1.0108689
# # lk.1_053 Starved    67786    0.8091367
# # lk.1_054     Fed    45454    0.9965984
# # lk.1_056     Fed    23812    1.0463683
# # lk.1_057     Fed    58775    1.0739499
# # lk.1_058     Fed    77134    1.0015224
# 
# # plotBCV after estimating common dispersion and tagwise dispersions
# #  in one run (recommended):
# design <- model.matrix(~group)
# y <- estimateDisp(y, design, robust=TRUE)
# pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/BCV.plot_12gLHLK_Train.Prot_20200706.pdf")
# # y <- estimateCommonDisp(y)
# # y <- estimateTrendedDisp(y)
# # y <- estimateTagwiseDisp(y)
# plotBCV(y)
# dev.off()
# 
# et <- exactTest(y, pair = c("Fed","Starved"))
# topTags(et, n = 20)
# # summary(decideTests(et)) # at 5% FDR
# # Comparison of groups:  Starved-Fed 
# # logFC    logCPM       PValue        FDR
# # CG5773      -2.231249  7.976927 5.935125e-06 0.03400826
# # fit         -2.256405  8.122176 3.069361e-05 0.08793720
# # CG33926     -2.912267  6.659562 1.188385e-04 0.22698162
# 
# # One gene (CG5773) coming out of default edgeR analysis from Starved - Fed comparison
# fit <- glmFit(y, design)
# lrt <- glmLRT(fit)
# # glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
# topTags(lrt, n = 30)
# summary(decideTests(lrt))
# # Coefficient:  groupStarved 
# # logFC    logCPM        LR       PValue        FDR
# # CG5773      -2.231249  7.976927 21.577092 3.398873e-06 0.01947554
# # fit         -2.256405  8.122176 17.290233 3.208316e-05 0.09191826
# # CG33926     -2.912267  6.659562 16.133157 5.904128e-05 0.11276884
# # CG1092      -2.946556  6.224946 11.979602 5.378607e-04 0.77048547
# # tim         -2.117463  7.110371 10.654482 1.098045e-03 0.98904635
# rownames(topTags(lrt, n = 30,p.value = 0.05)) -> edgeR.DE.gene
# # Three down-regulated genes (CG5773, fit, CG33926) coming out of the likelihood ratio test edgeR analysis from Starved - Fed comparison
# 
# 
# # DESeq2
# g.colData.df
# str(g.colData.df)
# as.factor(g.colData.df[,"Train_protocol"]) -> g.colData.df$Train_protocol
# dds <- DESeqDataSetFromMatrix(countData = g.count.df, 
#                               colData = g.colData.df, 
#                               design = ~ Train_protocol)
# # First, we filter out the lowly expressed genes, by removing those genes that 
# #   do not have at least 1 reads in at least half of the samples.
# dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
# filter <- rowSums(assay(dds)>1)>5
# table(filter) # FALSE 8860, TRUE 4170
# # Identify the 2000 most variable genes
# assay(dds.preF) %>% log1p %>% rowVars -> vars
# names(vars) <- rownames(dds.preF)
# vars <- sort(vars, decreasing = TRUE)
# head(vars)
# dds.preF <- dds.preF[names(vars)[1:2000],]
# dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
# dds.preF.analysis <- DESeq(dds.preF)
# res5.S2F <- results(dds.preF.analysis, alpha = 0.10, contrast = c("Train_protocol","Starved","Fed"))
# summary(res5.S2F)
# res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
# res5.S2F.Ordered
# 
# # log2 fold change (MLE): Train_protocol Starved vs Fed 
# # Wald test p-value: Train protocol Starved vs Fed 
# # DataFrame with 2000 rows and 6 columns
# # baseMean       log2FoldChange             lfcSE                 stat               pvalue               padj
# # <numeric>            <numeric>         <numeric>            <numeric>            <numeric>          <numeric>
# #   CG5773  10.6466395954202    -2.44152364404031 0.577392915385338    -4.22853065734431 2.35222467971256e-05 0.0470209713474541
# # CG33926 2.55716870427906    -3.20022223107304 0.938051360753999     -3.4115639771587 0.000645913417810402  0.645590461101497
# # Dh31    15.3490408422426    -3.35641179427909  1.03423093389026    -3.24532141158644   0.0011731812567957  0.652844938727196
# # fit     11.4100355125209    -2.35241301237106  0.73179422454291    -3.21458264287398  0.00130634304897888  0.652844938727196
# 
# # write.csv(as.data.frame(res5.S2F.Ordered),
# #           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk.1_6Starved-6Fed__20200706.csv")
# res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
# S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.050,])
# S2F.gene

#####
# # put all 18 LHLKs
# # DESeq2
# dds <- DESeqDataSetFromMatrix(countData = count.df[,-1], 
#                                 colData = colData.df[-1,], 
#                                 design = ~ Train_protocol)
# # First, we filter out the lowly expressed genes, by removing those genes that 
# #   do not have at least 1 reads in at least half of the samples.
# dds.preF <- dds[ rowSums(counts(dds)>1) >= ncol(counts(dds))/2,]
# filter <- rowSums(assay(dds)>1)>5
# table(filter) # FALSE 8600, TRUE 4430
# # Identify the 2000 most variable genes
# assay(dds.preF) %>% log1p %>% rowVars -> vars
# names(vars) <- rownames(dds.preF)
# vars <- sort(vars, decreasing = TRUE)
# head(vars)
# dds.preF <- dds.preF[names(vars)[1:2000],]
# dds.preF$Train_protocol <- relevel(dds.preF$Train_protocol, ref="Fed")
# dds.preF.analysis <- DESeq(dds.preF)
# res5.S2F <- results(dds.preF.analysis, alpha = 0.10, contrast = c("Train_protocol","Starved","Fed"))
# summary(res5.S2F)
# res5.S2F.Ordered <- res5.S2F[order(res5.S2F$padj),]
# res5.S2F.Ordered
# # write.csv(as.data.frame(res5.S2F.Ordered),
# #           file="/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/Keene_LHLK/analysis/lk.1_Starved-Fed_0.10_20200630.csv")
# res5.S2F.Ordered$padj <- ifelse(is.na(res5.S2F.Ordered$padj), 1, res5.S2F.Ordered$padj)
# S2F.gene = rownames(res5.S2F.Ordered[res5.S2F.Ordered$padj<=0.10,])
# S2F.gene


# Fed="#DFDF00", Starved="#009F7F"
df
gene.matrix["Hr38",]

-----

# relative overall heatmap
log2((1 + tpm.g.dpm.dir[,A.gDPM.id]) / (1 + apply(tpm.g.dpm.dir[,A.gDPM.id], 1, mean))) -> rel.tpm.A.g.dpm.dir
rel.tpm.A.g.dpm.dir
rel.tpm.A.g.dpm.dir[rel.tpm.A.g.dpm.dir < -2] <- -2
rel.tpm.A.g.dpm.dir[rel.tpm.A.g.dpm.dir > 2] <- 2
curBreaks <- seq(-2, 2, length.out=101)

pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/seq.data/DPM1-3.final/2019.rev/heatmap.d-cluster.noRownames.rel_dir.7A.DPMs.20190124.pdf", height = 10)
pheatmap(rel.tpm.A.g.dpm.dir, cluster_rows=TRUE, show_rownames=F, scale = "none",
         cluster_cols=F, annotation_col= NA, main = "Relative gene expression of 7 DPMs \n(log2 TPM/mean)",
         breaks = curBreaks, color = colorRampPalette(c("steelBlue2","white","darkOrange2"))(100),
         #         fontsize = 14, fontsize_col = 20,
         treeheight_row = 0, cellwidth = NA)
dev.off()

# get the MED detected gene count

mean(colData.g.dpm[A.gDPM.id,]$d.g.num)
colData.g.dpm[A.gDPM.id,]
# median of the detected gene num among the 31 good DPMs is 8205; among the 7 airflow DPMs is 7163
pdf("/Users/maxwellshih/Library/CloudStorage/Dropbox-JoshDubnau/Maxwell Shih/Maxwell_Lab/Memory/RNA-Seq_CEL-Seq/seq.data/DPM1-3.final/2019.rev/plot_bar_7A.DPMs_de.g.20200413.pdf")
ggplot(data = colData.g.dpm[A.gDPM.id,], aes(x = rownames(colData.g.dpm[A.gDPM.id,]), y=d.g.num)) + geom_bar(stat = "identity", width = 0.6) +
  scale_x_discrete(limits=rownames(colData.g.dpm[A.gDPM.id,]) ) + theme_gray(base_size = 20) +
  labs(title="", x ="Cell id", y = "Detected gene count") + geom_hline(yintercept=7163, linetype="dashed", color = "red") +
  geom_text(x = 5, y = 9000,label = paste('median = ',median(colData.g.dpm[A.gDPM.id,]$d.g.num), sep = ''), color="red", size = 10, family="Courier", fontface="plain")
dev.off()
#

