setwd("/home/zhuob/Project2014/Project1/data")



ds <- fit.nb.glm(nb.da, x)


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


data(arab)

pr <- 2 # average count per sample

# keep rows with at least mean 1 count 
arab1 <- arab1[rowSums(arab1)>pr*dim(arab1)[2], ]
arab2 <- arab2[rowSums(arab2)>pr*dim(arab2)[2], ]
arab6 <- arab6[rowSums(arab6)>pr*dim(arab6)[2], ]
arab7 <- arab7[rowSums(arab7)>pr*dim(arab7)[2], ]
arab8 <- arab8[rowSums(arab7)>pr*dim(arab7)[2], ]
arab9 <- arab9[rowSums(arab9)>pr*dim(arab9)[2], ]

arab0 <- arab[rowSums(arab)>pr*dim(arab)[2], ]


group1 <- c(1, 1, 1, 2, 2, 2)
group2 <- c(1, 1, 2, 2, 3, 3, 4, 4)
group6 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
group7 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
group8 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
group9 <- c(1, 1, 1, 2, 2, 2)
group0 <- c(1, 1, 1, 2, 2, 2)

head(nb.da$counts)

nf = estimate.norm.factors(arab1)
nb.da <- prepare.nb.data(arab1, norm.factors=nf)
nb.da$counts= as.matrix(nb.da$counts)
x = model.matrix(~ as.factor(group1))
dtt1 <- fit.nb.glm(nb.da, x)


nf = estimate.norm.factors(arab6)
nb.da <- prepare.nb.data(arab6, norm.factors=nf)
nb.da$counts= as.matrix(nb.da$counts)
x = model.matrix(~ as.factor(group6))
dtt6 <- fit.nb.glm(nb.da, x)


nf = estimate.norm.factors(arab7)
nb.da <- prepare.nb.data(arab7, norm.factors=nf)
nb.da$counts= as.matrix(nb.da$counts)
x = model.matrix(~ as.factor(group7))
dtt7 <- fit.nb.glm(nb.da, x)


nf = estimate.norm.factors(arab0)
nb.da <- prepare.nb.data(arab0, norm.factors=nf)
nb.da$counts= as.matrix(nb.da$counts)
x = model.matrix(~ as.factor(group0))
dtt0 <- fit.nb.glm(nb.da, x)


phi <- expandAsMatrix(dtt0$phi, c(dim(arab0)[1], dim(arab0)[2]))

system.time(sigma.0 <- estimate.dispersion.var.edgeR(nb.da, phi, x = x))




