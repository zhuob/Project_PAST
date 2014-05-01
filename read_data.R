

unwanted <- c("ambiguous","no_feature", "not_aligned", "too_low_aQual", "alignment_not_unique")

# 1-3 are control groups
arab.1 <- read.csv("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab1/arab1.csv", header=T)
arab1 <- arab.1[, 2:dim(arab.1)[2]]
row.names(arab1) <- arab1[,1]
arab1 <- arab1[, -1]
index <- which(row.names(arab1) %in% unwanted)
arab1 <- arab1[-index, ]
saveRDS(arab1, "/home/zhuob/Project2014/Project1/data/arab1.rds")


# 1-2 are control
arab.2 <- read.csv("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab2/arab2.csv", header=T)
arab2 <- arab.2[, 2:dim(arab.2)[2]]
index <- which(arab2$Gene %in% unwanted)
arab2[, 1]=substr(arab2[, 1], 1, 9)
row.names(arab2) <- arab2[,1]
arab2 <- arab2[-index, ]
arab2 <- arab2[(2:(dim(arab2)[1]-1)), -1]
## the "ambiguous" "no_feature" "not_aligned " "too_low_aQual" "alignment_not_unique"
saveRDS(arab2, "/home/zhuob/Project2014/Project1/data/arab2.rds")



## for arab6, if I check the website, then sample 5- 18 are used below
## Col-0 and Col-0 A are compared, and so C24 to C24 A
arab6 <- read.table("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab6/GSE42957_Read_Counts_TAIR10_genes.txt", header=T)
colnames(arab6)[1]="Gene"
index <- which(arab6$Gene %in% unwanted)
row.names(arab6) <- arab6[,1]
arab6 <- arab6[, -1]
saveRDS(arab6, "/home/zhuob/Project2014/Project1/data/arab6.rds")

# 10-12 are control group
arab7 <- read.table("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab7/GSE38879_sam2counts_RV.txt", header=T)
colnames(arab7)[1]="Gene"
index <- which(arab7$Gene %in% unwanted)
arab7[, 1]=substr(arab7[, 1], 1, 9)
row.names(arab7) <- arab7[,1]
arab7 <- arab7[, -1]
# delete NAs
arab7 <- arab7[complete.cases(arab7),] 
# head(arab7)
saveRDS(arab7, "/home/zhuob/Project2014/Project1/data/arab7.rds")

#  1-3 are control groups for arab8
#  the PLATFORM is Illumina Genome Analyzer IIx (Arabidopsis thaliana)
arab8 <- read.csv("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab8/arab8.csv", header=T)
arab8 <- arab8[, -1]
index <- which(arab8$Gene %in% unwanted)
row.names(arab8) <- arab8[,1]
arab8 <- arab8[, -1]
saveRDS(arab8, "/home/zhuob/Project2014/Project1/data/arab8.rds")


#  4-6 are control groups for arab9
arab.9 <- read.csv("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab9/arab9.csv", header=T)
arab9 <- arab.9[, 2:dim(arab.9)[2]]
index <- which(arab9$Gene %in% unwanted)
row.names(arab9) <- arab9[,1]
arab9 <- arab9[, -1]
saveRDS(arab9, "/home/zhuob/Project2014/Project1/data/arab9.rds")

arab10 <- read.table("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab10/GSE39463_CompleteCountTable_At.txt", header=T)
index <- which(row.names(arab1)%in% unwanted)
saveRDS(arab10, "/home/zhuob/Project2014/Project1/data/arab10.rds")



## column 1 is the gene name, col6 is the count data

setwd("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab12/")

## lots of other types, I screen type = mRNA

arab12.1 <- read.table("GSM595272_2_4cells_WT_unique_read_count.txt",header=T)[, c(1, 5,6)]
arab12.1 <- arab12.1[which(arab12.1[, 2]=="mRNA"),]
colnames(arab12.1)[3] = "WT_cells"

arab12.2 <- read.table("GSM595273_globular_WT_unique_read_count.txt",header=T)[, c(1,5, 6)]
arab12.2 <- arab12.2[which(arab12.2[, 2]=="mRNA"),]
colnames(arab12.2)[3] = "WT_globular"

arab12.3 <- read.table("GSM595274_2_4cells_kypKYP_unique_read_count.txt",header=T)[, c(1,5, 6)]
arab12.3 <- arab12.3[which(arab12.3[, 2]=="mRNA"),]
colnames(arab12.3)[3] = "kypKYP"

arab12.4 <- read.table("GSM838183_MC_LerCol_11_unique_read_count.txt",header=T)[, c(1,5, 6)]
arab12.4 <- arab12.4[which(arab12.4[, 2]=="mRNA"),]
colnames(arab12.4)[3] = "MC_LerCol_11"


arab12.5 <- read.table("GSM838184_MC_kypCol_12_unique_read_count.txt",header=T)[, c(1,5, 6)]
arab12.5 <- arab12.5[which(arab12.5[, 2]=="mRNA"),]
colnames(arab12.5)[3] = "MC_kypCol_12"



arab12.6 <- read.table("GSM949727_seedCoat_WT_unique_read_count.txt",header=T)[, c(1, 5,6)]
arab12.6 <- arab12.6[which(arab12.6[, 2]=="mRNA"),]
colnames(arab12.6)[3] = "seedCoat_WT"


merge1 <- merge(arab12.1[, c(1,3)], arab12.2[, c(1,3)], by="transcript")
merge2 <- merge(merge1, arab12.3[, c(1,3)], by="transcript")
merge3 <- merge(merge2, arab12.4[, c(1,3)], by="transcript")
merge4 <- merge(merge3, arab12.5[, c(1,3)], by="transcript")
arab12 <- merge(merge4, arab12.6[, c(1,3)], by="transcript")
colnames(arab12)[1]="Gene"
arab12[, 1]=substr(arab12[, 1], 1, 9)
arab12 <- arab12[which(!duplicated(arab12$Gene)), ]
index <- which(arab8$Gene %in% unwanted)
row.names(arab12) <- arab12[,1]
arab12 <- arab12[, -1]
saveRDS(arab12, "/home/zhuob/Project2014/Project1/data/arab12.rds")



