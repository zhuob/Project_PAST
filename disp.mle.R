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
dim(arab2)
length(obj2$Gene)
obj1 <- pi.disp.matrix(arab1, group1, tol.phi = 1e-7) # if I don't filter out, convergence issue
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

obj <- join_all(list(obj1, obj2, obj6, obj7, obj8, obj9, obj0), by= "Gene")
colnames(obj) <- c("Gene", "dis1", "pi1", "disp6", "pi6", "disp7", "pi7",
                   "disp8", "pi8", "disp9", "pi9", "disp0", "pi0")
colnames(obj) <- c("Gene", "dis1", "pi1", "disp2", "pi2", "disp6", "pi6", "disp7", 
                   "pi7","disp8", "pi8", "disp9", "pi9", "disp0", "pi0")

obj <- obj[complete.cases(obj), ]
head(obj)
xx1  <- log(obj[, 3])
xx2  <- log(obj[, 5]) 
xx6  <- log(obj[, 7])
xx7  <- log(obj[, 9])
xx8  <- log(obj[, 11])
xx9  <- log(obj[, 13])
xx0  <- log(obj[, 15])


yy1 <- obj[, 2]
yy2 <- obj[, 4]
yy6 <- obj[, 6]
yy7 <- obj[, 8]
yy8 <- obj[, 10]
yy9 <- obj[, 12]
yy0 <- obj[, 14]

## the following has convergence issue
# deg = 1
# fit <- glm(yy0 ~ poly(xx0, degree = deg) + poly(xx1, degree=deg) +
#              poly(xx6, degree=deg) + poly(xx7, degree=deg) +
#              poly(xx8, degree=deg) + poly(xx9, degree=deg) 
#            , family=Gamma(link="log"))
deg =2


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



nf <- estimate.norm.factors(lh1$counts, method="AH2010")
nb.da <- prepare.nb.data(lh1$counts)
nb.da$counts <- as.matrix(nb.da$counts)
x <- model.matrix(~as.factor(group0))
disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ")


sigma.0 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x) #1.000719
sigma0 <- sigma.est(obj=lh1) #  0.5588521
sigma1 <- sigma.est(obj=lh2) #  0.5418044
disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ", method="ML")
sigma.ML.0 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x)

phi <- list()
phi$estimates <- expandAsMatrix(yy0, c(dim(obj)[1], length(group0))) 
sigma.gen.0 <- estimate.dispersion.var(nb.da, dispersion = phi, x = x)


##  arab1 ##--------------------------------------------
fit5 <- glm(yy1~poly(xx1, degree=deg)+poly(xx6, degree=deg)
             , family=Gamma(link="log"))
summary(fit5)

fit6 <- glm(yy1~poly(xx1, degree=deg)+poly(xx6, degree=deg)+ log(yy6)
            , family=Gamma(link="log"))


lh3 <- likelihood.score(data= arab1, phi=fitted(fit5), dat1=obj, group=group1)
lh4 <- likelihood.score(data= arab1, phi=fitted(fit6), dat1=obj, group=group1)

nf <- estimate.norm.factors(lh3$counts, method="AH2010")
nb.da <- prepare.nb.data(lh3$counts)
nb.da$counts <- as.matrix(nb.da$counts)
x <- model.matrix(~as.factor(lh3$group))
disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ")


sigma.1 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x) # 1.746364
sigma2 <- sigma.est(obj=lh3) #  0.9700425 
sigma3 <- sigma.est(obj=lh4) #  0.7277991


disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ", method="ML")
sigma.ML.1 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x)

phi <- list()
phi$estimates <- expandAsMatrix(yy1, c(dim(obj)[1], length(group1))) 
sigma.gen.1 <- estimate.dispersion.var(nb.da, dispersion = phi, x = x)




##  arab6 ####
fit9 <- glm(yy6~poly(xx0, degree= deg)+ poly(xx1, degree=deg)+poly(xx6, degree=deg)
            + poly(xx7, degree=deg)+ poly(xx8, degree=deg)+ poly(xx9, degree=deg), 
            family=Gamma(link="log"))

fit10 <- glm(yy6~poly(xx0, degree= deg)+ poly(xx1, degree=deg)+poly(xx6, degree=deg)
             + poly(xx7, degree=deg)+ poly(xx8, degree=deg)+ poly(xx9, degree=deg)
             + log (yy7), 
             family=Gamma(link="log"))


lh10$likelihood
lh9 <- likelihood.score(data= arab6, phi=fitted(fit9), dat1=obj, group=group6)
lh10 <- likelihood.score(data= arab6, phi=fitted(fit10), dat1=obj, group=group6)

lh <- lh9
nf <- estimate.norm.factors(lh$counts, method="AH2010")
nb.da <- prepare.nb.data(lh$counts)
nb.da$counts <- as.matrix(nb.da$counts)
x <- model.matrix(~as.factor(lh$group))
disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ")


sigma.3 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x) # 2.697607
sigma6 <- sigma.est(obj=lh9) # 1.113342
sigma7 <- sigma.est(obj=lh10) # 0.9961775

disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ", method="ML")
sigma.ML.6 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x)

phi <- list()
phi$estimates <- expandAsMatrix(yy6, c(dim(obj)[1], length(group6))) 
sigma.gen.6 <- estimate.dispersion.var(nb.da, dispersion = phi, x = x)



# arab7
fit7 <- glm(yy7~poly(xx0, degree= deg)+ poly(xx1, degree=deg)+poly(xx6, degree=deg)
            + poly(xx7, degree=deg)+ poly(xx8, degree=deg)+ poly(xx9, degree=deg), 
            family=Gamma(link="log"))



fit8 <- glm(yy7~poly(xx0, degree= deg)+ poly(xx1, degree=deg)+poly(xx6, degree=deg)
            + poly(xx7, degree=deg)+ poly(xx8, degree=deg)+ poly(xx9, degree=deg)
             + log (yy1), 
            family=Gamma(link="log"))


lh6$likelihood
lh5 <- likelihood.score(data= arab7, phi=fitted(fit7), dat1=obj, group=group7)
lh6 <- likelihood.score(data= arab7, phi=fitted(fit8), dat1=obj, group=group7)


nf <- estimate.norm.factors(lh6$counts, method="AH2010")
nb.da <- prepare.nb.data(lh6$counts)
nb.da$counts <- as.matrix(nb.da$counts)
x <- model.matrix(~as.factor(lh6$group))
disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ")


sigma.2 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x) # 1.115855
sigma4 <- sigma.est(obj=lh5) #  0.4454039
sigma5 <- sigma.est(obj=lh6) #  0.3917444


disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ", method="ML")
sigma.ML.7 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x)

phi <- list()
phi$estimates <- expandAsMatrix(yy7, c(dim(obj)[1], length(group7))) 
sigma.gen.7 <- estimate.dispersion.var(nb.da, dispersion = phi, x = x)



## arab8 ## ----
fit11 <- glm(yy8~poly(xx0, degree= deg)+ poly(xx1, degree=deg)+poly(xx6, degree=deg)
            + poly(xx7, degree=deg)+ poly(xx8, degree=deg)+ poly(xx9, degree=deg), 
            family=Gamma(link="log"))

fit12 <- glm(yy8~poly(xx0, degree= deg)+ poly(xx1, degree=deg)+poly(xx6, degree=deg)
             + poly(xx7, degree=deg)+ poly(xx8, degree=deg)+ poly(xx9, degree=deg)
             + log (yy7), 
             family=Gamma(link="log"))

fit13 <- glm(yy8~poly(xx0, degree= deg)+ poly(xx1, degree=deg)+poly(xx6, degree=deg)
             + poly(xx7, degree=deg)+ poly(xx8, degree=deg)+ poly(xx9, degree=deg)
             + log (yy1), 
             family=Gamma(link="log"))



lh13$likelihood
lh11 <- likelihood.score(data= arab8, phi=fitted(fit11), dat1=obj, group=group8)
lh12 <- likelihood.score(data= arab8, phi=fitted(fit12), dat1=obj, group=group8)

lh <- lh11
nf <- estimate.norm.factors(lh$counts, method="AH2010")
nb.da <- prepare.nb.data(lh$counts)
nb.da$counts <- as.matrix(nb.da$counts)
x <- model.matrix(~as.factor(lh$group))
disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ")


sigma.4 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x) # 1.689247
sigma8 <- sigma.est(obj=lh11) # 1.04766
sigma9 <- sigma.est(obj=lh12) # 1.042108




disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ", method="ML")
sigma.ML.8 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x)

phi <- list()
phi$estimates <- expandAsMatrix(yy8, c(dim(obj)[1], length(group8))) 
sigma.gen.8 <- estimate.dispersion.var(nb.da, dispersion = phi, x = x)


lh13 <- likelihood.score(data= arab8, phi=fitted(fit13), dat1=obj, group=group8)
sigma10 <- sigma.est(obj=lh13) # 1.035786


round(cor(log(obj[, c(2, 4, 6, 8, 10, 12)])),3)

setwd("/home/zhuob/Project2014/Project1/data")

mag <- readRDS("mag.rds")
id <- apply(mag[,-1]>1e-9, 1, all)
mag.n <- mag[id,]
dim(mag.n) # 1801 15


## ##    estimating dispersion by NBQ and obtain sigma^2


nf <- estimate.norm.factors(arab6, method="AH2010")
nb.da <- prepare.nb.data(arab6)
nb.da$counts <- as.matrix(nb.da$counts)
x <- model.matrix(~as.factor(group6))
disp.arab0 <- estimate.dispersion(nb.da, x=x, model="NBQ", method= "MAPL")


re.6 <- estimate.dispersion.var(nb.da, dispersion=disp.arab0, x = x) # 1.689247
