
setwd("/home/zhuob/Project2014/Project1/data")
library(edgeR)
library(NBPSeq)
library(SeqDisp)
library(dplyr)


arab1 <- readRDS("arab1.rds")
arab6 <- readRDS("arab6.rds")


group.lab6 = c(1,1, 2,2,3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
mod.po <- poly.gamma.reg(counts= arab6, x = group.lab6)

group.id <- c(1, 1, 1, 2, 2, 2)
mod2 <- poly.gamma.reg(counts = arab1, x = group.id)


source("/home/zhuob/Project2014/Project1/R/dispersion.edgeR.R")
source("/home/zhuob/Project2014/Project1/R/create.geneColumn.R")

disp6 <- dispersion.edgeR(arab6, group.lab6)
disp1 <- dispersion.edgeR(arab1, group.id)

y <- prepare.nb.data(arab6)
rel.freq6 <- data.frame(mean.rel.freq6=apply(y$rel.frequencies, 1, mean))


y <- prepare.nb.data(arab1)
rel.freq1 <- data.frame(mean.rel.freq1=apply(y$rel.frequencies, 1, mean))


head(rel.freq6)
disp6 <- create.geneColumn(disp6)
disp1 <- create.geneColumn(disp1)
rel.freq6 <- create.geneColumn(rel.freq6)
rel.freq1 <- create.geneColumn(rel.freq1)

target <- join_all(list(disp6, rel.freq6, disp1, rel.freq1), by= "Gene")
names(target)[-1]= c("disp6","rel.freq6", "disp1","rel.freq1")
dat <- target[complete.cases(target),] # remove columns containing NAs


id <- which( (dat[,5]>0) & (dat[, 3]>0) )
id2 <- which(all(dat[,-1]>0))

dat <- dat[id, ]
rho <- round( cor(log(dat[, -1])), 3) ## correlation coefficient



NBPSeq:::smart.plot(log(dat[,c(5, 2)]), clip=16, pch =20, main =paste("cor = ", rho[1,4]))
NBPSeq:::smart.plot(log(dat[,c(3, 2)]), clip=16, pch =20, main =paste("cor = ", rho[1,2]))
NBPSeq:::smart.plot(log(dat[,c(4, 2)]), clip=16, pch =20, main =paste("cor = ", rho[1,3]))








x1 <- log(dat$disp1)
x2 <- log(dat$rel.freq6)
x3 <- log(dat$rel.freq1)
y <- dat$disp6

#####     pairs(log(dat[,-1]), pch=20)

mod1 <- glm(y ~x1+ x2+ x3 , family=Gamma(link="log"))
mod2 <- glm(y~x1 + x2, family=Gamma(link="log"))
mod3 <- glm(y~ x2, family=Gamma(link="log"))
mod4 <- glm(y~1, family=Gamma(link="log"))

### quadratic term 
mod5 <- glm(y ~x1 + poly(x2, degree =2)+ x3 , family=Gamma(link="log"))
mod6 <- glm(y ~poly(x1, degree=2) + poly(x2, degree =2)+ x3 , family=Gamma(link="log"))
mod7 <- glm(y ~poly(x1, degree=2) + poly(x2, degree =2)+ poly(x3, degree=2)
            , family=Gamma(link="log"))

### cubic term
mod8 <- glm(y ~poly(x1, degree=2) + poly(x2, degree =3)+ poly(x3, degree=2)
            , family=Gamma(link="log"))

mod9 <- glm(y ~poly(x1, degree=3) + poly(x2, degree =3)+ poly(x3, degree=2)
            , family=Gamma(link="log"))

mod10 <- glm(y ~poly(x1, degree=2) + poly(x2, degree =3)+ poly(x3, degree=3)
            , family=Gamma(link="log"))

mod11 <- glm(y ~poly(x1, degree=3) + poly(x2, degree =3)+ poly(x3, degree=3)
             , family=Gamma(link="log"))

summary(mod10)

mod12 <- glm(y ~ poly(x2, degree =3)+ poly(x3, degree=3)
             , family=Gamma(link="log"))

## newly added 
mod13 <- glm(dat$disp1 ~ poly(x2, degree =3)+ poly(x3, degree=3)
             , family=Gamma(link="log"))


summary(mod11)
AIC(mod10, mod12)

AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8,mod9, mod10, mod11)


## plot(mod10)

plot(x1, resid(mod10))
NBPSeq:::smart.plot(resid(mod12),x1, clip=32,pch= 20)
cor(resid(mod12), x1)

## see if there is pattern
NBPSeq:::smart.plot(resid(mod13),resid(mod12), clip=32,pch= 20)
abline(0,1)




