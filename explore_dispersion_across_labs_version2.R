
setwd("/home/zhuob/Project2014/Project1/data")
library(edgeR)
library(NBPSeq)
library(dplyr)

source("/home/zhuob/Project2014/Project1/R/disp.mom.R")

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

source("/home/zhuob/Project2014/Project1/R/dispersion.edgeR.R")
source("/home/zhuob/Project2014/Project1/R/create.geneColumn.R")


group.id <- c(1, 1, 1, 2, 2, 2)
disp.lab1.cox_reid <- dispersion.edgeR(arab1, group.id)
head(disp.lab1.cox_reid)
disp.lab1.cox_reid <- create.geneColumn(disp.lab1.cox_reid)

group <- c(1,1,2,2,3,3, 4, 4)
disp.lab2.cox_reid <- dispersion.edgeR(arab2, group)
disp.lab2.cox_reid <- create.geneColumn(disp.lab2.cox_reid)


group.id = c(1,1, 2,2,3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
disp.lab6.cox_reid <- dispersion.edgeR(arab6[, 1:14], group.id)
disp.lab6.cox_reid <- create.geneColumn(disp.lab6.cox_reid)



## data <- merge(disp.lab1.cox_reid, disp.lab6.cox_reid, by ="row.names")
## cor(log(data[, 2:3]))


group = c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4)
disp.lab7.cox_reid <- dispersion.edgeR(arab7, group)
disp.lab7.cox_reid <- create.geneColumn(disp.lab7.cox_reid)


group = c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4)
disp.lab8.cox_reid <- dispersion.edgeR(arab8, group)
disp.lab8.cox_reid <- create.geneColumn(disp.lab8.cox_reid)


group = c(1,1, 1,2,2, 2)
disp.lab9.cox_reid <- dispersion.edgeR(arab9, group)
disp.lab9.cox_reid <- create.geneColumn(disp.lab9.cox_reid)


group = c(1,1, 2, 2, 3,3)
disp.lab12.cox_reid <- dispersion.edgeR(arab12, group)
disp.lab12.cox_reid <- create.geneColumn(disp.lab12.cox_reid)


disp.7labs.cox_reid <- join_all(
    list(disp.lab1.cox_reid, disp.lab2.cox_reid,
        disp.lab6.cox_reid, disp.lab7.cox_reid,disp.lab8.cox_reid,
        disp.lab9.cox_reid, disp.lab12.cox_reid),
      by = "Gene")


colnames(disp.7labs.cox_reid) <- c("Gene", "lab1", "lab2", "lab6",
                                  "lab7", "lab8", "lab9", "lab12")


disp.7labs.cox_reid <- disp.7labs.cox_reid[complete.cases(disp.7labs.cox_reid),]

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

colSums(arab1[, -1])
colSums(arab2[, -1])
colSums(arab6[, -1])
colSums(arab7[, -1])
colSums(arab8[, -1])
colSums(arab9[, -1])
colSums(arab12[, -1])

rank(c(
     mean(colSums(arab1[, -1])),
     mean(colSums(arab2[, -1])),
     mean(colSums(arab6[, -1])),
     mean(colSums(arab7[, -1])),
     mean( colSums(arab8[, -1])),
     mean(colSums(arab9[, -1])),
     mean( colSums(arab12[, -1]))  
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






