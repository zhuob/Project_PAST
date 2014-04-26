

# 1-3 are control groups
arab.1 <- read.csv("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab1/arab1.csv", header=T)
arab1 <- arab.1[, 2:dim(arab.1)[2]]

# 1-2 are control
arab.2 <- read.csv("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab2/arab2.csv", header=T)
arab2 <- arab.2[, 2:dim(arab.2)[2]]
arab2[, 1]=substr(arab2[, 1], 1, 9)

arab6 <- read.table("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab6/GSE42957_Read_Counts_TAIR10_genes.txt", header=T)
colnames(arab6)[1]="Gene"

# 10-12 are control group
arab7 <- read.table("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab7/GSE38879_sam2counts_RV.txt", header=T)
colnames(arab7)[1]="Gene"
arab7[, 1]=substr(arab7[, 1], 1, 9)

#  1-3 are control groups for arab8
#  the PLATFORM is Illumina Genome Analyzer IIx (Arabidopsis thaliana)
arab8 <- read.csv("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab8/arab8.csv", header=T)
#  4-6 are control groups for arab9
arab.9 <- read.csv("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab9/arab9.csv", header=T)
arab9 <- arab.9[, 2:dim(arab.9)[2]]
arab10 <- read.table("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab10/GSE39463_CompleteCountTable_At.txt", header=T)

## column 1 is the gene name, col6 is the count data

setwd("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab12/")

arab12.1 <- read.table("GSM595272_2_4cells_WT_unique_read_count.txt",header=T)[, c(1, 6)]
colnames(arab12.1)[2] = "WT1"
arab12.2 <- read.table("GSM595273_globular_WT_unique_read_count.txt",header=T)[, c(1, 6)]
colnames(arab12.2)[2] = "WT2"
arab12.3 <- read.table("GSM595274_2_4cells_kypKYP_unique_read_count.txt",header=T)[, c(1, 6)]
colnames(arab12.3)[2] = "MC_LerCol_11"

arab12.4 <- read.table("GSM838183_MC_LerCol_11_unique_read_count.txt",header=T)[, c(1, 6)]
colnames(arab12.4)[2] = "MC_LerCol_12"

arab12.5 <- read.table("GSM838184_MC_kypCol_12_unique_read_count.txt",header=T)[, c(1, 6)]
colnames(arab12.5)[2] = "MC_kypCol_12"

arab12.6 <- read.table("GSM949727_seedCoat_WT_unique_read_count.txt",header=T)[, c(1, 6)]
colnames(arab12.6)[2] = "seedCoat_WT"
merge1 <- merge(arab12.1, arab12.2, by="transcript")
merge2 <- merge(merge1, arab12.3, by="transcript")
merge3 <- merge(merge2, arab12.4, by="transcript")
merge4 <- merge(merge3, arab12.5, by="transcript")
arab12 <- merge(merge4, arab12.6, by="transcript")
colnames(arab12)[1]="Gene"
arab12[, 1]=substr(arab12[, 1], 1, 9)


setwd("/home/zhuob/Project2014/arabidopsis_data")
write.table(arab1, "arab1.txt")
write.table(arab2, "arab2.txt")
write.table(arab6, "arab6.txt")
write.table(arab7, "arab7.txt")
write.table(arab8, "arab8.txt")
write.table(arab9, "arab9.txt")
write.table(arab10, "arab10.txt")
write.table(arab12, "arab12.txt")


# setwd("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab/arab13/")
# arab13.1 <- read.table("GSM607723_Col-dmexLer-WT_Endosperm_RNA-Seq_cDNA_mapping_imprinting.txt", header=T)
# head(arab13.1)


# # merge data example ## 
# df1 = data.frame(CustomerId=c(1:6),Product=c(rep("Toaster",3),rep("Radio",3)))
# df2 = data.frame(CustomerId=c(2,4,6),State=c(rep("Alabama",2),rep("Ohio",1)))
# 
# df1
# df2
# merge(x = df1, y = df2, by = "CustomerId", all = TRUE) # outer join
# merge(x = df1, y = df2, by = "CustomerId", all.x=TRUE) # left join
# merge(x = df1, y = df2, by = "CustomerId", all.y=TRUE) # right outer
# merge(x = df1, y = df2, by = NULL)  # cross join

arab.1<- merge(x = arab1[, 2:5], y = arab2[,2:4], by="Gene" )

arab.no.NA<-  arab.1[complete.cases(arab.1),]
arab.no.NA[1:20, ]

arab.2 <- merge(x=arab.1, y=arab8[, 2:5], by="Gene")
arab.3 <- merge(x=arab.2, y=arab9[, c(2, 6:8)], by="Gene")
head(arab.3)
library(NBPSeq)
data(arab)


write.csv(arab.3, "C:/Users/zhuob/Dropbox/Zhuo/FindData/arab/arabCtrl.csv")

# combine the arabidopsis data in NBPSeq
#  name the unnamed first column

s1 <- apply(arab.3[,2:4],1, mean)
s2 <- apply(arab.3[, 5:6], 1, mean)
s3 <- apply(arab.3[, 7:9], 1, mean)
s4 <- apply(arab.3[, 10:12], 1, mean)

pdf("C:/Users/zhuob/Dropbox/Zhuo/FindData/arab/log_mean.pdf",width=7,height=5)
par(mfrow=c(2, 3))
plot(log(s1), log(s2))
plot(log(s1), log(s3))
plot(log(s1), log(s4))
plot(log(s2), log(s3))
plot(log(s2), log(s4))
plot(log(s3), log(s4))
dev.off()



head(arab9)

arab.4 <- merge(x=arab.2, y=arab9[, c(2, 3:8)], by="Gene")

dat <- arab.4[, 2:15]
dat.1 <- dat + 1
boxplot(log(dat.1))


pairs(log(arab9[, 3:8]+1))












