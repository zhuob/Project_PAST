
setwd("/home/zhuob/Project2014/Project1/data")
library(edgeR)
library(NBPSeq)
library(SeqDisp)
library(dplyr)

source("/home/zhuob/Project2014/Project1/R/dispersion.edgeR.R")
source("/home/zhuob/Project2014/Project1/R/create.geneColumn.R")
source("/home/zhuob/Project2014/Project1/R/mean.disp.matrix.R")
source("/home/zhuob/Project2014/Project1/R/likelihood.score.R")
source("/home/zhuob/Project2014/Project1/R/estimate.dispersion.var.R")

arab1 <- readRDS("arab1.rds")
arab2 <- readRDS("arab2.rds")
arab6 <- readRDS("arab6.rds")
arab7 <- readRDS("arab7.rds")
arab8 <- readRDS("arab8.rds")
arab9 <- readRDS("arab9.rds")

head(arab2)
pr <- 5 # average count per sample

# keep rows with at least mean 1 count 
arab1 <- arab1[rowSums(arab1)>pr*dim(arab1)[2], ]
arab2 <- arab2[rowSums(arab2)>pr*dim(arab2)[2], ]
arab6 <- arab6[rowSums(arab6)>pr*dim(arab6)[2], ]
arab7 <- arab7[rowSums(arab7)>pr*dim(arab7)[2], ]
arab8 <- arab8[rowSums(arab7)>pr*dim(arab7)[2], ]
arab9 <- arab9[rowSums(arab9)>pr*dim(arab9)[2], ]

group1 <- c(1, 1, 1, 2, 2, 2)
group2 <- c(1, 1, 2, 2, 3, 3, 4, 4)
group6 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
group7 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
group8 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
group9 <- c(1, 1, 1, 2, 2, 2)



dtt1 <- mean.disp.matrix(arab1, group1)
dtt2 <- mean.disp.matrix(arab2, group2)
dtt6 <- mean.disp.matrix(arab6, group6)
dtt7 <- mean.disp.matrix(arab7, group7)
dtt8 <- mean.disp.matrix(arab8, group8)
dtt9 <- mean.disp.matrix(arab9, group9)


dat1 <- join_all(list(dtt1, dtt2, dtt6, dtt7, dtt8, dtt9), by= "Gene")
colnames(dat1) <- c("Gene", "disp1","pi1", "disp2","pi2", "disp6", "pi6", 
                   "disp7", "pi7", "disp8","pi8", "disp9","pi9")

dat1 <- dat1[complete.cases(dat1), ]
dim(dat1)



round(cor(log(dat1[, -1])), 2)
# it seems group 2 has consistently lower correlation with other groups

# # again if we sample, we might reach different conclusions
# set.seed(20)
# n <- 1000
# id <- sample(1:dim(dat)[1], n)
# dat <- dat[id, ]


xx1 <- log(dat1$pi1)
xx2 <- log(dat1$pi2)
xx6 <- log(dat1$pi6)
xx7 <- log(dat1$pi7)
xx8 <- log(dat1$pi8)
xx9 <- log(dat1$pi9)


yy1 <- dat1$disp1
yy2 <- dat1$disp2
yy6 <- dat1$disp6
yy7 <- dat1$disp7
yy8 <- dat1$disp8
yy9 <- dat1$disp9

library(MASS)

# y1, x1***, x2***, x6***, x7.    x8**   x9***
# y2,        x2***,        x7*  , x8***
# y6, x1***,        x6***, x7***, x8***, x9***
# y7, x1***, x2***, x6***, x7***, x8***, x9***
# y8, x1***,        x6*,          x8***, x9***
# y9, x1***, x2**, x6***, x7***



## function to calculate likelihood


# without dispersion term

str(arab2)

undebug(likelihood.score)

l0 <- likelihood.score(data= arab1, pred="NULL",response= yy1, d=1, group=group1) #  -597254.7
l1 <- likelihood.score(data= arab1, pred=yy6, response= yy1, d=1, group=group1) # -595482


source("/home/zhuob/Project2014/Project1/R/sigma.est.R")

s1.1 <- sigma.est(l1, method="NBQ")
s1.2 <- sigma.est(l1, method="edgeR")
s1.3 <- sigma.est(l1, method= "gamma.reg")


l6 <- likelihood.score(data= arab6, pred=yy7, response= yy6, d=1, group=group6)
s6.1 <- sigma.est(l6, method="NBQ")
s6.2 <- sigma.est(l6, method="edgeR")
s6.3 <- sigma.est(l6, method= "gamma.reg")


l7 <- likelihood.score(data= arab1, pred=yy7, response= yy1, d=1, group=group1) 
s7.1 <- sigma.est(l7, method="NBQ")
s7.2 <- sigma.est(l7, method="edgeR")
s7.3 <- sigma.est(l7, method= "gamma.reg")



# figure out the deviance 


# quantify.noise(dt="mouse", m=5000, sub.col=1:6, evaluate=c(TRUE, FALSE, FALSE), seed=539, path.o=path.mus)

# dispersion = estimate.dispersion(nb.data = nb.data, x = x, model = "NBQ", 
# method = "MAPL")
# head(dispersion$estimates)
# [,1]        [,2]        [,3]        [,4]        [,5]        [,6]
# [1,] 0.061066466 0.061066466 0.061066466 0.055306168 0.055306168 0.055306168
# [2,] 0.003503290 0.003503290 0.003503290 0.003485469 0.003485469 0.003485469