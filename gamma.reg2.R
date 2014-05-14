
setwd("/home/zhuob/Project2014/Project1/data")
library(edgeR)
library(NBPSeq)
library(SeqDisp)
library(dplyr)

source("/home/zhuob/Project2014/Project1/R/dispersion.edgeR.R")
source("/home/zhuob/Project2014/Project1/R/create.geneColumn.R")
source("/home/zhuob/Project2014/Project1/R/mean.disp.matrix.R")

arab1 <- readRDS("arab1.rds")
arab6 <- readRDS("arab6.rds")
arab7 <- readRDS("arab7.rds")

# keep rows with at least mean 1 count 
arab1 <- arab1[rowSums(arab1)>dim(arab1)[2], ]
arab6 <- arab6[rowSums(arab6)>dim(arab6)[2], ]
arab7 <- arab7[rowSums(arab7)>dim(arab7)[2], ]


group1 <- c(1, 1, 1, 2, 2, 2)
group6 <- c(1,1, 2,2,3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
group7 <- c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4)



dtt1 <- mean.disp.matrix(arab1, group1)
dtt6 <- mean.disp.matrix(arab6, group6)
dtt7 <- mean.disp.matrix(arab7, group7)

## use join_all function in dplyr package

dat <- join_all(list(dtt1, dtt6, dtt7), by= "Gene")
colnames(dat) <- c("Gene", "disp1","pi1", "disp6", "pi6", "disp7", "pi7")
dim(dat)
index <- which(dat=="NaN")
tail(dat)
dat <- dat[complete.cases(dat), ]
dim(dat)

x1 <- log(dat$pi1)
x6 <- log(dat$pi6)
x7 <- log(dat$pi7)

y1 <- dat$disp1
y6 <- dat$disp6
y7 <- dat$disp7


m1 <- glm(y1 ~x1+ x6 + x7 , family=Gamma(link="log"))
n1 <- glm(y6 ~x1+ x6 + x7 , family=Gamma(link="log"))
p1 <- glm(y7 ~x1+ x6 + x7 , family=Gamma(link="log"))

# m7 has the smallesst AIC with all terms significant
m4 <- glm(y1 ~poly(x1, degree= 2)+ x6 + poly(x7, degree= 2) , family=Gamma(link="log"))
m5 <- glm(y1 ~poly(x1, degree= 2)+ x6 + poly(x7, degree= 1) , family=Gamma(link="log"))
m6 <- glm(y1 ~poly(x1, degree= 3)+ x6 + poly(x7, degree= 1) , family=Gamma(link="log"))
m7 <- glm(y1 ~poly(x1, degree= 3) + poly(x7, degree= 1) , family=Gamma(link="log"))
m8 <- glm(y1 ~poly(x1, degree= 3) + poly(x7, degree= 2) , family=Gamma(link="log"))


AIC(m1, m4, m5, m6, m7, m8)

summary(m7)

n2 <- glm(y6 ~poly(x1, degree= 3) + poly(x6, degree= 3) + 
            poly(x7, degree= 2) , family=Gamma(link="log"))
n3 <- glm(y6 ~poly(x1, degree= 3) + poly(x6, degree= 3) + 
            poly(x7, degree= 1) , family=Gamma(link="log"))
n4 <- glm(y6 ~poly(x1, degree= 3) + poly(x6, degree= 3) 
           , family=Gamma(link="log"))

AIC(n1, n2, n3, n4)
summary(n4)


p2 <- glm(y7 ~poly(x1, degree= 3) + poly(x6, degree= 3) + 
            poly(x7, degree= 3) , family=Gamma(link="log"))
p3 <- glm(y7 ~poly(x1, degree= 2) + poly(x6, degree= 3) + 
            poly(x7, degree= 3) , family=Gamma(link="log"))
p4 <- glm(y7 ~poly(x1, degree= 1) + poly(x6, degree= 3) + 
            poly(x7, degree= 3) , family=Gamma(link="log"))

AIC(p1, p2, p3, p4)
summary(p4)



par(mfrow=c(1, 1))
NBPSeq:::smart.plot(resid(p4),resid(n4), clip=16,pch= 20,
                    main="lab2 vs lab3 (y = lab2)")
abline(0,1)
NBPSeq:::smart.plot(resid(m7),resid(n4), clip=16,pch= 20,
                    main="lab1 vs lab2")
abline(0,1)
NBPSeq:::smart.plot(resid(p4),resid(m7), clip=16,pch= 20,
                    main="lab1 vs lab3")
abline(0,1)
# mtext("dispersion regressed over rel.freq", side = 3, line = -20, outer = TRUE)
# dev.off()


# p <- ggplot(residual, aes(resid1,resid2))
# p + geom_point(alpha= 1/40)

# ggplot(residual, aes(resid1,resid2, color=resid1)) +
# geom_point() + scale_colour_gradientn(colours = rainbow(5)) + theme_bw()


# ggplot(residual, aes(resid1,resid2))+ geom_density2d()


## get a subsample to see if the result is consistent

## -----  don't know if this is a problem that we concern?

set.seed(20)
n <- 1000
id <- sample(1:dim(dat)[1], n)
dat1 <- dat[id, ]

x1 <- log(dat1$pi1)
x6 <- log(dat1$pi6)
x7 <- log(dat1$pi7)

y1 <- dat1$disp1
y6 <- dat1$disp6
y7 <- dat1$disp7

p1 <- glm(y7 ~x1+ x6 + x7 , family=Gamma(link="log"))
p2 <- glm(y7 ~poly(x1, degree= 3) + poly(x6, degree= 3) + 
            poly(x7, degree= 3) , family=Gamma(link="log"))
p3 <- glm(y7 ~poly(x1, degree= 2) + poly(x6, degree= 3) + 
            poly(x7, degree= 3) , family=Gamma(link="log"))
p4 <- glm(y7 ~poly(x1, degree= 1) + poly(x6, degree= 3) + 
            poly(x7, degree= 3) , family=Gamma(link="log"))

AIC(p1, p2, p3, p4)
summary(p2)

n2 <- glm(y6 ~poly(x1, degree= 3) + poly(x6, degree= 3) + 
            poly(x7, degree= 2) , family=Gamma(link="log"))
n3 <- glm(y6 ~poly(x1, degree= 3) + poly(x6, degree= 3) + 
            poly(x7, degree= 1) , family=Gamma(link="log"))
n4 <- glm(y6 ~poly(x1, degree= 3) + poly(x6, degree= 3) 
          , family=Gamma(link="log"))




