setwd("/home/zhuob/Project2014/arabidopsis_data/")
 
#  read a bulk of data simultaneously 
# series <- c(1, 2, 6, 7, 8, 9, 10 ,12)
# datalist= list()
# 
# for ( i in 1: length(series))
# {
#   datalist[[i]] <- read.table(paste("arab", series[i],".txt", sep=""), header=T)
#   
# }

arab1 <- read.table("arab1.txt", header = T) 
## make sure at least one count on average for each row
arab1 <- arab1[ which(rowSums(arab1[, 2:dim(arab1)[2]])>= dim(arab1)[2]),] 

arab2 <- read.table("arab2.txt", header = T)
arab2 <- arab2[ which(rowSums(arab2[, 2:dim(arab2)[2]])>= dim(arab2)[2]),]

arab6 <- read.table("arab6.txt", header = T)
arab6 <- arab6[ which(rowSums(arab6[, 2:dim(arab6)[2]])>= dim(arab6)[2]),]

arab7 <- read.table("arab7.txt", header = T)
arab7 <- arab7[ which(rowSums(arab7[, 2:dim(arab7)[2]])>= dim(arab7)[2]),]

arab8 <- read.table("arab8.txt", header = T)
arab8 <- arab8[ which(rowSums(arab8[, 2:dim(arab8)[2]])>= dim(arab8)[2]),]

arab9 <- read.table("arab9.txt", header = T)
arab9 <- arab9[ which(rowSums(arab9[, 2:dim(arab9)[2]])>= dim(arab9)[2]),]

arab10 <- read.table("arab10.txt", header = T)
arab10 <- arab10[ which(rowSums(arab10[, 2:dim(arab10)[2]])>= dim(arab10)[2]),]

arab12 <- read.table("arab12.txt", header = T)
arab12 <- arab12[ which(rowSums(arab12[, 2:dim(arab12)[2]])>= dim(arab12)[2]),]

## Method of moments:  var / mean^2
arab1$mean = apply(as.matrix(arab1[,2:dim(arab1)[2]]), 1,  mean)
arab1$var = apply(as.matrix(arab1[,2:dim(arab1)[2]]), 1,  var)
arab1$var.devideby.squaredmean.lab1 = arab1$var/arab1$mean^2
arab1.disp <- arab1[, c(1, dim(arab1)[2])]

arab2$mean = apply(as.matrix(arab2[,2:dim(arab2)[2]]), 1,  mean)
arab2$var = apply(as.matrix(arab2[,2:dim(arab2)[2]]), 1,  var)
arab2$var.devideby.squaredmean.lab2 = arab2$var/arab2$mean^2
arab2.disp <- arab2[, c(1, dim(arab2)[2])]

arab6$mean = apply(as.matrix(arab6[,2:dim(arab6)[2]]), 1,  mean)
arab6$var = apply(as.matrix(arab6[,2:dim(arab6)[2]]), 1,  var)
arab6$var.devideby.squaredmean.lab6 = arab6$var/arab6$mean^2
arab6.disp <- arab6[, c(1, dim(arab6)[2])]


arab7$mean = apply(as.matrix(arab7[,2:dim(arab7)[2]]), 1,  mean)
arab7$var = apply(as.matrix(arab7[,2:dim(arab7)[2]]), 1,  var)
arab7$var.devideby.squaredmean.lab7 = arab7$var/arab7$mean^2
arab7.disp <- arab7[, c(1, dim(arab7)[2])]

arab8$mean = apply(as.matrix(arab8[,2:dim(arab8)[2]]), 1,  mean)
arab8$var = apply(as.matrix(arab8[,2:dim(arab8)[2]]), 1,  var)
arab8$var.devideby.squaredmean.lab8 = arab8$var/arab8$mean^2
arab8.disp <- arab8[, c(1, dim(arab8)[2])]


arab9$mean = apply(as.matrix(arab9[,2:dim(arab9)[2]]), 1,  mean)
arab9$var = apply(as.matrix(arab9[,2:dim(arab9)[2]]), 1,  var)
arab9$var.devideby.squaredmean.lab9 = arab9$var/arab9$mean^2
arab9.disp <- arab9[, c(1, dim(arab9)[2])]

arab12$mean = apply(as.matrix(arab12[,2:dim(arab12)[2]]), 1,  mean)
arab12$var = apply(as.matrix(arab12[,2:dim(arab12)[2]]), 1,  var)
arab12$var.devideby.squaredmean.lab12 = arab12$var/arab12$mean^2
arab12.disp <- arab12[, c(1, dim(arab12)[2])]

merge1 <- merge(arab1.disp, arab2.disp, by = "Gene")
merge2 <- merge(merge1, arab6.disp, by = "Gene")
merge3 <- merge(merge2, arab7.disp, by = "Gene")
merge4 <- merge(merge3, arab8.disp, by = "Gene")
merge5 <- merge(merge4, arab9.disp, by = "Gene")
disp.7labs <- merge(merge5, arab12.disp, by = "Gene")



### using edgeR estimateTagwiseDisp() function qCML method to obtain dispersion



#----------------------------------------------------------

library(edgeR)
arab1.edgeR <- arab1[, 2:7]
group = factor(c(1,1,1,2,2,2))
y <- DGEList(counts=arab1.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateTagwiseDisp(y)
disp.lab1 <- y$tagwise.dispersion
disp.lab1.edgeR <- data.frame(arab1$Gene, disp.lab1)
colnames(disp.lab1.edgeR)[1]= "Gene"

arab2.edgeR <- arab2[, 2:9]
group = factor(c(1,1,2,2,3,3, 4, 4))
y <- DGEList(counts=arab2.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateTagwiseDisp(y)
disp.lab2 <- y$tagwise.dispersion
disp.lab2.edgeR <- data.frame(arab2$Gene, disp.lab2)
colnames(disp.lab2.edgeR)[1]= "Gene"

names(arab6)
## group factors see http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42957
arab6.edgeR <- arab6[, 2:15]
group = factor (c(1,1, 2,2,3, 3, 4, 4, 5, 5, 6, 6, 7, 7))
y <- DGEList(counts=arab6.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateTagwiseDisp(y)
disp.lab6 <- y$tagwise.dispersion
disp.lab6.edgeR <- data.frame(arab6$Gene, disp.lab6)
colnames(disp.lab6.edgeR)[1]= "Gene"

arab7.edgeR <- arab7[, 2:13]
group = factor (c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4))
y <- DGEList(counts=arab7.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateTagwiseDisp(y)
disp.lab7 <- y$tagwise.dispersion
disp.lab7.edgeR <- data.frame(arab7$Gene, disp.lab7)
colnames(disp.lab7.edgeR)[1]= "Gene"


arab8.edgeR <- arab8[, 2:13]
group = factor (c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4))
y <- DGEList(counts=arab8.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateTagwiseDisp(y)
disp.lab8 <- y$tagwise.dispersion
disp.lab8.edgeR <- data.frame(arab8$Gene, disp.lab8)
head(disp.lab8.edgeR)
colnames(disp.lab8.edgeR)[1]= "Gene"

arab9.edgeR <- arab9[, 2:7]
group = factor (c(1,1, 1,2,2, 2))
y <- DGEList(counts=arab9.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateTagwiseDisp(y)
disp.lab9 <- y$tagwise.dispersion
disp.lab9.edgeR <- data.frame(arab9$Gene, disp.lab9)
colnames(disp.lab9.edgeR)[1]= "Gene"

arab12.edgeR <- arab12[, 2:7]
group = factor (c(1,1, 1, 2, 2,2))  # ti should actually be 1 rep per group. !!
y <- DGEList(counts=arab12.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateTagwiseDisp(y)
disp.lab12 <- y$tagwise.dispersion
disp.lab12.edgeR <- data.frame(arab12$Gene, disp.lab12)
colnames(disp.lab12.edgeR)[1]= "Gene"

merge1 <- merge(disp.lab1.edgeR, disp.lab2.edgeR, by = "Gene")
merge2 <- merge(merge1, disp.lab6.edgeR, by = "Gene")
merge3 <- merge(merge2, disp.lab7.edgeR, by = "Gene")
merge4 <- merge(merge3, disp.lab8.edgeR, by = "Gene")
merge5 <- merge(merge4, disp.lab9.edgeR, by = "Gene")
#  dispersion matrix estimated by qCML in edgeR
disp.7labs.qCML <- merge(merge5, disp.lab12.edgeR, by = "Gene")





# -- -------------------------  Cox-Reid Profile Likelihood ------------------------

arab1.edgeR <- arab1[, 2:7]
group <- factor(c(1,1,1,2,2,2))
design <- model.matrix(~group)
y <- DGEList(counts=arab1.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
y <- estimateGLMTagwiseDisp(y, design)
disp.lab1.cr <- y$tagwise.dispersion
disp.lab1.cox_reid <- data.frame(arab1$Gene, disp.lab1.cr)
colnames(disp.lab1.cox_reid)[1]= "Gene"



arab2.edgeR <- arab2[, 2:9]
group = factor(c(1,1,2,2,3,3, 4, 4))
design <- model.matrix(~group)
y <- DGEList(counts=arab2.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
y <- estimateGLMTagwiseDisp(y, design)
disp.lab2.cr <- y$tagwise.dispersion
disp.lab2.cox_reid <- data.frame(arab2$Gene, disp.lab2.cr)
colnames(disp.lab2.cox_reid)[1]= "Gene"


## group factors see http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42957
arab6.edgeR <- arab6[, 2:15]
group = factor (c(1,1, 2,2,3, 3, 4, 4, 5, 5, 6, 6, 7, 7))
y <- DGEList(counts=arab6.edgeR,group=group)
design <- model.matrix(~group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
y <- estimateGLMTagwiseDisp(y, design)

disp.lab6.cr <- y$tagwise.dispersion
disp.lab6.cox_reid <- data.frame(arab6$Gene, disp.lab6.cr)
colnames(disp.lab6.cox_reid)[1]= "Gene"

arab7.edgeR <- arab7[, 2:13]
group = factor (c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4))
design <- model.matrix(~group)
y <- DGEList(counts=arab7.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
y <- estimateGLMTagwiseDisp(y, design)
disp.lab7.cr <- y$tagwise.dispersion
disp.lab7.cox_reid <- data.frame(arab7$Gene, disp.lab7.cr)
colnames(disp.lab7.cox_reid)[1]= "Gene"


arab8.edgeR <- arab8[, 2:13]
group = factor (c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4))
design <- model.matrix(~group)
y <- DGEList(counts=arab8.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
y <- estimateGLMTagwiseDisp(y, design)
disp.lab8.cr <- y$tagwise.dispersion
disp.lab8.cox_reid <- data.frame(arab8$Gene, disp.lab8.cr)
head(disp.lab8.cr)
colnames(disp.lab8.cox_reid)[1]= "Gene"

arab9.edgeR <- arab9[, 2:7]
group = factor (c(1,1, 1,2,2, 2))
design = model.matrix(~group)
y <- DGEList(counts=arab9.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
y <- estimateGLMTagwiseDisp(y, design)
disp.lab9.cr <- y$tagwise.dispersion
disp.lab9.cox_reid <- data.frame(arab9$Gene, disp.lab9.cr)
colnames(disp.lab9.cox_reid)[1]= "Gene"

arab12.edgeR <- arab12[, 2:7]
group = factor (c(1,1, 2, 2, 3,3))  # ti should actually be 1 rep per group. !!
design <- model.matrix(~group)
y <- DGEList(counts=arab12.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
y <- estimateGLMTagwiseDisp(y, design)
disp.lab12.cr <- y$tagwise.dispersion
disp.lab12.cox_reid <- data.frame(arab12$Gene, disp.lab12.cr)
colnames(disp.lab12.cox_reid)[1]= "Gene"

merge1 <- merge(disp.lab1.cox_reid, disp.lab2.cox_reid, by = "Gene")
merge2 <- merge(merge1, disp.lab6.cox_reid, by = "Gene")
merge3 <- merge(merge2, disp.lab7.cox_reid, by = "Gene")
merge4 <- merge(merge3, disp.lab8.cox_reid, by = "Gene")
merge5 <- merge(merge4, disp.lab9.cox_reid, by = "Gene")
#  dispersion matrix estimated by qCML in edgeR
disp.7labs.cox_reid <- merge(merge5, disp.lab12.cox_reid, by = "Gene")





dim(disp.7labs.qCML)
dim(disp.7labs)
dim(disp.7labs.cox_reid)



##################
##  plot the results

lab <- c(1, 2, 6, 7, 8, 9, 12)

n <- 100  ## number of genes to be displayed

set.seed(120)
index <- sample(1:dim(disp.7labs.qCML)[1],n)

## create any label I want in x-axis
xval <- seq(1:7)
plot(xval, log( disp.7labs[1,-1]), xaxt="n", xlab= "lab",
     ylab= "var/mean^2", pch = 20, ylim= c(-5,3),
     main=paste("# of genes displayed:", n))
# change x-axis with specified lab labels
axis(1, at=1:7, labels=lab)

lines(xval, log(disp.7labs[1,-1]), col= 6)
for ( i in 1: n)
{
  points(xval, log(disp.7labs[index[i], -1]), pch =20)
  lines(xval, log(disp.7labs[index[i], -1]), col= i) 
}

###
colMeans(disp.7labs[, 2:dim(disp.7labs)[2]])



plot(xval, log( disp.7labs.qCML[1,-1]), xaxt="n", xlab= "lab",
     ylab= "dispersion in log scale", pch = 20, ylim= c(-8,1),
     main=paste("# of genes displayed:", n, "qCML(edgeR)"))
# change x-axis with specified lab labels
axis(1, at=1:7, labels=lab)

lines(xval, log(disp.7labs.qCML[1,-1]), col= 6)
for ( i in 1: n)
{
  points(xval, log(disp.7labs.qCML[index[i], -1]), pch =20)
  lines(xval, log(disp.7labs.qCML[index[i], -1]), col= i) 
}

###
colMeans(disp.7labs.qCML[, 2:dim(disp.7labs.qCML)[2]])

#-------------------Cox-Reid  profile likelihood ----------

plot(xval, log( disp.7labs.cox_reid[1,-1]), xaxt="n", xlab= "lab",
     ylab= "dispersion in log scale", pch = 20, ylim= c(-8,1),
     main=paste("# of genes displayed:", n, "(cox-reid (edgeR))"))
# change x-axis with specified lab labels
axis(1, at=1:7, labels=lab)

lines(xval, log(disp.7labs.cox_reid[1,-1]), col= 6)
for ( i in 1: n)
{
  points(xval, log(disp.7labs.cox_reid[index[i], -1]), pch =20)
  lines(xval, log(disp.7labs.cox_reid[index[i], -1]), col= i) 
}

###
colMeans(disp.7labs.cox_reid[, 2:dim(disp.7labs.cox_reid)[2]])



colSums(arab1[, 2:7])
colSums(arab2[, 2:9])
colSums(arab6[, 2:13])
colSums(arab7[, 2:13])
colSums(arab8[, 2:13])
colSums(arab9[, 2:7])
colSums(arab12[, 2:7])

rank(c(
     mean(colSums(arab1[, 2:7])),
     mean(colSums(arab2[, 2:9])),
     mean(colSums(arab6[, 2:15])),
     mean(colSums(arab7[, 2:13])),
     mean( colSums(arab8[, 2:13])),
     mean(colSums(arab9[, 2:7])),
     mean( colSums(arab12[, 2:7]))  
)
  )

round(cor(log(disp.7labs.cox_reid[, -1])),2)
pairs(log(disp.7labs.cox_reid[, -1]))
head(disp.7labs.cox_reid)
 
NBPSeq:::smart.plot(log(disp.7labs.cox_reid[,c(4,2)]), clip = 16, pch = 20
                    , xlab="log-dispersion (lab1)", ylab="log-dispersion (lab 6)",
        main= paste( "cor =", round(cor(log(disp.7labs.cox_reid[, c(2, 4)]))[2, 1], 2)))

head(arab6)
arab6$rowsum.lab6 = rowSums(arab6[,1:12])
arab6$Gene <- row.names(arab6)
arab6.rowsum.disp = merge(disp.lab7.cox_reid, arab6[, c(dim(arab6)[2]-1, dim(arab6)[2])], by= "Gene")
head(arab6[, c(1, dim(arab6)[2]-1)])
head(arab6.rowsum.disp)


head(disp.6labs.cox_reid)

round(cor(log(disp.7labs.cox_reid[, c(2, 4)]))[2, 1], 2)
round(cor(log(arab6.rowsum.disp[, 2:3]))[2, 1], 2)
NBPSeq:::smart.plot(log(arab6.rowsum.disp[, 2:3]), clip= 16, pch = 20, 
  main= paste( "both in log scale, cor= ",
               round(cor(log(arab6.rowsum.disp[, 2:3]))[2, 1], 2)))



