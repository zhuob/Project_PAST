
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


dtt1 <- mean.disp.matrix(arab1, group1) #0.001
dtt2 <- mean.disp.matrix(arab2, group2) #0.0001
dtt6 <- mean.disp.matrix(arab6, group6) #0.1
dtt7 <- mean.disp.matrix(arab7, group7) #0.001
dtt8 <- mean.disp.matrix(arab8, group8) #0.01
dtt9 <- mean.disp.matrix(arab9, group9) # 0.01
dtt0 <- mean.disp.matrix(arab0, group0) #0.1


dat1 <- join_all(list(dtt1, dtt2, dtt6, dtt7, dtt8, dtt9, dtt0), by= "Gene")
colnames(dat1) <- c("Gene", "disp1","pi1", "disp2","pi2", "disp6", "pi6", 
                   "disp7", "pi7", "disp8","pi8", "disp9","pi9", "disp0", "pi0")

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
xx0 <- log(dat1$pi0)

yy1 <- dat1$disp1
yy2 <- dat1$disp2
yy6 <- dat1$disp6
yy7 <- dat1$disp7
yy8 <- dat1$disp8
yy9 <- dat1$disp9
yy0 = dat1$disp0

library(MASS)

# y1, x1***, x2***, x6***, x7.    x8**   x9***
# y2,        x2***,        x7*  , x8***
# y6, x1***,        x6***, x7***, x8***, x9***
# y7, x1***, x2***, x6***, x7***, x8***, x9***
# y8, x1***,        x6*,          x8***, x9***
# y9, x1***, x2**, x6***, x7***





quantify.noise(dt="arab0", m=5000, sub.col=1:6,
            evaluate=c(TRUE, FALSE, FALSE), seed=539, path.o=path.mus)
quantify.noise(dt="arab0", m=5000, sub.col=1:6, 
            evaluate=c(FALSE, TRUE, FALSE), seed=539, path.o=path.mus)




fit <- glm(yy0 ~ poly(xx0, degree=deg), family=Gamma(link="log"))
phi.hat <- fitted(fit)
phi <- expandAsMatrix(phi.hat, c(dim(dat1)[1], dim(dat1)[2]))
length(phi.hat)


## function to calculate likelihood
source("/home/zhuob/Project2014/Project1/R/likelihood.score.R")



####  use arabidopsis data 

ll0 <- likelihood.score(data= arab0, pred="NULL",response= yy0, d=2, group=group1)

s.0 = sigma.est(ll0, method="gamma.reg") # 56.71848

ll1 <- likelihood.score(data= arab0, pred=yy1 ,response= yy0, d=2, group=group1)
s.1 = sigma.est(ll1, method="gamma.reg") # 2.482933
s.2 <- sigma.est(ll1, method="NBQ") #  1.457224
s.3 <- sigma.est(ll1, method="edgeR") #  0.7488082



l0 <- likelihood.score(data= arab1, pred="NULL",response= yy1, d=1, group=group1) #  -597254.7
l1 <- likelihood.score(data= arab1, pred=yy6, response= yy1, d=1, group=group1) # -595482

source("/home/zhuob/Project2014/Project1/R/sigma.est.R")

s1.0 = sigma.est(l0, method="gamma.reg") # 3.87
s1.1 <- sigma.est(l1, method="NBQ") # 2.87
s1.2 <- sigma.est(l1, method="edgeR") # 1.17
s1.3 <- sigma.est(l1, method= "gamma.reg") # 3.61



l6 <- likelihood.score(data= arab6, pred=yy1, response= yy6, d=1, group=group6)
s6.1 <- sigma.est(l6, method="NBQ") # 4.15
s6.2 <- sigma.est(l6, method="edgeR") # 2.68
s6.3 <- sigma.est(l6, method= "gamma.reg") # 6.34


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