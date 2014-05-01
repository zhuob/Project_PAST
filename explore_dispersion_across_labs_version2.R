
setwd("/home/zhuob/Project2014/Project1/data")
library(edgeR)
library(NBPSeq)

source("/home/zhuob/Project2014/Project1/man/disp.mom.R")

arab1 <- readRDS("arab1.rds")

arab2 <- readRDS("arab2.rds")
arab6 <- readRDS("arab6.rds")
arab7 <- readRDS("arab7.rds")
arab8 <- readRDS("arab8.rds")
arab9 <- readRDS("arab9.rds")
arab12 <- readRDS("arab12.rds")


arab1.disp <- dispersion.mom(counts= arab1, g= 2, r =3 )
arab2.disp <- dispersion.mom(arab2, g=4, r=2)
arab6.disp <- dispersion.mom(arab6, g=6, r=2)
arab7.disp <- dispersion.mom(arab7, g=4, r=3)
arab8.disp <- dispersion.mom(arab8, g=4, r=3)
arab9.disp <- dispersion.mom(arab9, g=2, r=3)
arab12.disp <- dispersion.mom(arab12, g=3, r=2)

merge <- merge(arab1.disp, arab6.disp, by = "row.names")
cor(log(merge[, -1]))

# 
# merge1 <- merge(arab1.disp, arab2.disp, by = "row.names")
# colnames(merge1)[1] = "Gene"
# head(merge1)
# 
# merge2 <- merge(merge1, arab6.disp, by = "Gene")
# merge3 <- merge(merge2, arab7.disp, by = "row.names")
# merge4 <- merge(merge3, arab8.disp, by = "row.names")
# merge5 <- merge(merge4, arab9.disp, by = "row.names")
# disp.7labs <- merge(merge5, arab12.disp, by = "row.names")


# -- -------------------------  Cox-Reid Profile Likelihood ------------------------

source("/home/zhuob/Project2014/Project1/man/disp.edgeR.R")


group.id <- c(1, 1, 1, 2, 2, 2)
disp.lab1.cox_reid <- dispersion.edgeR(arab1, group.id)

group <- c(1,1,2,2,3,3, 4, 4)
disp.lab2.cox_reid <- dispersion.edgeR(arab2, group)

group.id = c(1,1, 2,2,3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
disp.lab6.cox_reid <- dispersion.edgeR(arab6[, 1:14], group.id)



## data <- merge(disp.lab1.cox_reid, disp.lab6.cox_reid, by ="row.names")
## cor(log(data[, 2:3]))


group = c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4)
disp.lab7.cox_reid <- dispersion.edgeR(arab7, group)

group = c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4)
disp.lab8.cox_reid <- dispersion.edgeR(arab8, group)

group = c(1,1, 1,2,2, 2)
disp.lab9.cox_reid <- dispersion.edgeR(arab9, group)

group = c(1,1, 2, 2, 3,3)
disp.lab12.cox_reid <- dispersion.edgeR(arab12, group)


merge1 <- merge(disp.lab1.cox_reid, disp.lab2.cox_reid, by = "row.names")


## notice after Merge, the row.names become a new column
colnames(merge1)=c("Gene","disp.lab1", "disp.lab2")
disp.lab6.cox_reid[, 2]= row.names(disp.lab6.cox_reid) 
colnames(disp.lab6.cox_reid)=c( "disp.lab6","Gene")
merge2 <- merge(merge1, disp.lab6.cox_reid, by = "Gene")


disp.lab7.cox_reid[, 2]= row.names(disp.lab7.cox_reid) 
colnames(disp.lab7.cox_reid)=c( "disp.lab7","Gene")
merge3 <- merge(merge2, disp.lab7.cox_reid, by = "Gene")

disp.lab8.cox_reid[, 2]= row.names(disp.lab8.cox_reid) 
colnames(disp.lab8.cox_reid)=c( "disp.lab8","Gene")
merge4 <- merge(merge3, disp.lab8.cox_reid, by = "Gene")

disp.lab9.cox_reid[, 2]= row.names(disp.lab9.cox_reid) 
colnames(disp.lab9.cox_reid)=c( "disp.lab9","Gene")
merge5 <- merge(merge4, disp.lab9.cox_reid, by = "Gene")


disp.lab12.cox_reid[, 2]= row.names(disp.lab12.cox_reid) 
colnames(disp.lab12.cox_reid)=c( "disp.lab12","Gene")
disp.7labs.cox_reid <- merge(merge5, disp.lab12.cox_reid, by = "Gene")


head(disp.7labs.cox_reid)
dim(disp.7labs.cox_reid)


round(cor(log(disp.7labs.cox_reid[, -1])), 3)


##################
##  plot the results

lab <- c(1, 2, 6, 7, 8, 9, 12)

n <- 100  ## number of genes to be displayed

set.seed(120)
index <- sample(1:dim(disp.7labs.cox_reid)[1],n)

## create any label I want in x-axis
xval <- seq(1:length(lab))

#-------------------Cox-Reid  profile likelihood ----------

plot(xval, log( disp.7labs.cox_reid[1,-1]), xaxt="n", xlab= "lab",
     ylab= "dispersion in log scale", pch = 20, ylim= c(-8,1),
     main=paste("# of genes displayed:", n, "(cox-reid (edgeR))"))
# change x-axis with specified lab labels
axis(1, at=1:length(lab), labels=lab)

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

round(cor(log(disp.7labs.cox_reid[, -1])),3)
## pairs(log(disp.6labs.cox_reid[, -1]))


rho <- round(cor(log(disp.7labs.cox_reid[, c(2, 4)]))[2, 1], 3)
NBPSeq:::smart.plot(log(disp.7labs.cox_reid[,c(4,2)]), clip = 16, pch = 20
                    , xlab="log-dispersion (lab1)", ylab="log-dispersion (lab 6)",
        main= paste( "cor =", rho))


y <- prepare.nb.data(arab6)
rel.freq <-  apply(y$rel.frequencies, 1, sum)
## keep the gene names
rel.freq <- data.frame(rel.freq)

head(disp.lab6.cox_reid)
head(rel.freq)
y <- merge(disp.lab6.cox_reid, rel.freq, by ="row.names")
# remove rows with 0 rel. freq
head(y)
arab6.rel.freq.disp <- y[which(y[, 4]>0),-3] 

rho2 <-  round(cor(log(arab6.rel.freq.disp[, 2:3]))[2, 1], 3)
NBPSeq:::smart.plot(log(arab6.rel.freq.disp[, 2:3]), clip= 16, pch = 20, 
  main= paste( "both in log scale, cor= ",rho2))







library(SeqDisp)

group = c(1,1, 2,2,3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
x <- model.matrix(~group)
a <- poly.gamma.reg(counts= arab6, group)


counts <- arab6
head(counts)
m = dim(counts)[1]
n = dim(counts)[2]
mu.mom = expandAsMatrix(rowMeans(counts), dim=c(m,n))



