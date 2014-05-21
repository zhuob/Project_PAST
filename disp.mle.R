library(edgeR)
library(NBPSeq)
library(SeqDisp)
setwd("/home/zhuob/Project2014/Project1/data")


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
group1 <- c(1, 1, 1, 2, 2, 2)
group2 <- c(1, 1, 2, 2, 3, 3, 4, 4)
group6 <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
group7 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
group8 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
group9 <- c(1, 1, 1, 2, 2, 2)
group0 <- c(1, 1, 1, 2, 2, 2)


obj1 <- pi.disp.matrix(arab1, group1, tol.phi= 1e-7)
obj2 <- pi.disp.matrix(arab2, group2, tol.phi= 1e-7)
obj6 <- pi.disp.matrix(arab6, group6, tol.phi= 1e-7)
obj7 <- pi.disp.matrix(arab7, group7, tol.phi= 1e-7)
obj8 <- pi.disp.matrix(arab8, group8, tol.phi= 1e-7)
obj9 <- pi.disp.matrix(arab9, group9, tol.phi= 1e-7)
obj0 <- pi.disp.matrix(arab0, group0, tol.phi= 1e-7)



# ma0 <- data.frame(rf0 = rel.freq0, disp0 = dtt0$phi)


# library(edgeR)
# phi <- expandAsMatrix(dtt0$phi, c(dim(arab0)[1], dim(arab0)[2]))

# system.time(sigma.0 <- estimate.dispersion.var(nb.da, dtt0, x = x))

obj <- join_all(list(obj1, obj6, obj7, obj8, obj9, obj0), by= "Gene")
obj <- obj[complete.cases(obj), ]
dim(obj)

xx1  <- log(obj[, 3])
xx6  <- log(obj[, 5])
xx7  <- log(obj[, 7])
xx8  <- log(obj[, 9])
xx9  <- log(obj[, 11])
xx0  <- log(obj[, 13])


yy1 <- obj[, 2]
yy6 <- obj[, 4]
yy7 <- obj[, 6]
yy8 <- obj[, 8]
yy9 <- obj[, 10]
yy0 <- obj[, 12]

## the following has convergence issue
# deg = 1
# fit <- glm(yy0 ~ poly(xx0, degree = deg) + poly(xx1, degree=deg) +
#              poly(xx6, degree=deg) + poly(xx7, degree=deg) +
#              poly(xx8, degree=deg) + poly(xx9, degree=deg) 
#            , family=Gamma(link="log"))

fit2 <- glm(yy0~poly(xx0, degree= deg)+ poly(xx1, degree=deg)+poly(xx6, degree=deg)
            + poly(xx7, degree=deg)+ poly(xx8, degree=deg), 
            family=Gamma(link="log"))

deg =2
fit3 <- glm(yy0~poly(xx0, degree= deg)+ poly(xx7, degree=deg), 
            family=Gamma(link="log"))

summary(fit3)

fit4 <- glm(yy0~poly(xx0, degree= deg)+ poly(xx7, degree=deg) + log(yy7), 
            family=Gamma(link="log"))

summary(fit4)

lh1 <- likelihood.score(data= arab0, phi=fitted(fit3), dat1=obj, group=group0)
lh2 <- likelihood.score(data= arab0, phi=fitted(fit4), dat1=obj, group=group0)

sigma0 <- sigma.est(obj=lh1)
sigma1 <- sigma.est(obi=lh2)







setwd("/home/zhuob/Project2014/Project1/data")

mag <- readRDS("mag.rds")
id <- apply(mag[,-1]>1e-9, 1, all)
mag.n <- mag[id,]
dim(mag.n) # 1801 15
# 
# xx1 <- log(mag$rf1)
# xx2 <- log(mag$rf2)
# xx6 <- log(mag$rf6)
# xx7 <- log(mag$rf7)
# xx8 <- log(mag$rf8)
# xx9 <- log(mag$rf9)
# xx0 <- log(mag$rf0)
# 
# yy1 <- mag$disp1
# yy2 <- mag$disp2
# yy6 <- mag$disp6
# yy7 <- mag$disp7
# yy8 <- mag$disp8
# yy9 <- mag$disp9
# yy0 = mag$disp0
# 
# NBPSeq:::smart.plot(log(mag[, 5]), log(mag[, 9]), clip=12, xlab ="arab1", ylab="arab6" )
# 
# par(mfrow=c(2,2))
# hist(log(mag[, 3]), main ="arab")
# hist(log(mag[,5]), main="arab1")
# hist(log(mag[,7]), main="arab2")
# hist(log(mag[,9]), main="arab6")
# #dev.off()
# log(colMeans(mag[, -1]))
# 
# 
# hist(log(mag[, 11]), main ="arab7")
# hist(log(mag[,13]), main="arab8")
# hist(log(mag[,15]), main="arab9")
# 
# group0
# library(MASS)
# # debug(likelihood.score)
# # undebug(glm)
# # undebug(glm.fit)
# ll0 <- likelihood.score(data= arab1, pred=yy6, response= yy1, d=2, group=group1)
# # undebug(likelihood.score)
# 
# str(y)
# str(xx0)
# deg = 2
#  y = yy1[1:100]
# fit <- glm(y ~ poly(xx0, degree = deg) + poly(xx1, degree=deg) +
#              poly(xx2, degree=deg) +
#              poly(xx6, degree=deg) + poly(xx7, degree=deg) +
#              poly(xx8, degree=deg) + poly(xx9, degree=deg)
#            , family=Gamma(link="log"))




