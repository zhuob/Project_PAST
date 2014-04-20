
## set the working directory
setwd("/home/zhuob/Dropbox/Zhuo/Research/Project_2014/Data/arab")

arab <- read.csv("arabCtrl.csv", header=T)
head(arab)
arab.dat <- arab[, 3:dim(arab)[2]] 

## package glmmADMB deals with count data via glmm
library(glmmADMB)

### get a random subset of arab.dat to run this glmm model
#   treat gene id as a fixed effect, and group as random effect.

group <- as.factor(c(1, 1,1, 2,2,3,3,3, 4, 4, 4))

## filter genes with 1 count per sample on average
arab.filter <- as.matrix(arab.dat[rowSums(arab.dat)>= dim(arab.dat)[2], ])

#  number of genes to be sampled
n <- 100
set.seed(128) 

set.seed(5)

## generate the sequence of genes to be sampled
id <- sample (1:dim(arab.filter)[1], n)  

arab.sample <- arab.filter[id, ]

##  store the variances for gene sampled 
vari <- c()


for (i in 1:n){
  y <- arab.sample[i, ]
  a <- glmmadmb(y ~1, random=  ~1 | group, zeroInflation=F, 
                link="log",  family = "nbinom")
  vari[i] <- as.numeric(a$S)  #### output the variance
}

###  some warning message would pop up. 
#  1: In glmmadmb(y ~ 1, random = ~1 | group, zeroInflation = F, link = "log",  :
# Convergence failed:log-likelihood of gradient= -0.00148958
# convergence issue 

# try y <- c(34, 43, 42, 7, 3, 0, 0, 0, 0, 0, 0 )
# or y <- arab.filter[21251,] 


## draw a histogram of the log variances.
hist(log(vari), main=paste("number of genes =", n))

id2 <- which(log(vari)< -10)
## to see which genes have estimated variance of 0
arab.pre[id[id2],]


# we can have glmmPQL estimate the variance, but it cannot be output.
y <- arab.pre[id[25],]
b <- glmmPQL(y ~1, random=  ~1 | group, family = quasipoisson)

a <- glmmadmb(y~1, random=  ~1 | group, zeroInflation=F, 
              link="log",  family = "nbinom")
summary(a)
## in general, a gives smaller variance than b







## resluts for the following genes are compared to NLMIXED in SAS and they are pretty close

# head(arab.dat)
# y1 <- as.numeric(arab.dat[6, ])
# y2 <- as.numeric(arab.dat[7502, ])
# y3 <- as.numeric(arab.dat[19895, ])
# y4 <- as.numeric(arab.pre[id[25],])
# y5 <- as.numeric(arab.dat[23,])
# y6 <- as.numeric(arab.dat[id[12],])
# y7 <- as.numeric(arab.dat[7782,])
# y8 <- as.numeric(arab.dat[12475,])
# cbind(y1, y2, y3, y4, y5, y6, y7, y8)

# a <- glmmadmb(y1~1, random=  ~1 | group, zeroInflation=F, 
#               link="log",  family = "nbinom")

#  installing R-INLA package
# source("http://www.math.ntnu.no/inla/givemeINLA-testing.R") 

# to check how glmmadmb  approximate the likelihood
file.show(system.file("tpl","glmmadmb.tpl",package="glmmADMB"))




##  to see whether it makes difference if library sizes are normalized
library(NBPSeq)


norm.factors <- estimate.norm.factors(arab.pre)

lib.size <- colSums(arab.pre)
max(norm.factors)/min(norm.factors) #  3.22
max(lib.size)/min(lib.size) # 6.44


## initializing relative frequency matrixa
arab.filter.norm<- matrix(0, ncol= dim(arab.pre)[2], nrow = dim(arab.pre)[1])


## multiply each count by the corresponding norm.factor and round to the 
## nearest integer
for ( i in 1: dim(arab.filter)[2])
{
  arab.filter.norm[, i]= round(arab.filter[, i]*norm.factors[i], 0)
  
}

arab.sample.norm <- arab.filter.norm[id, ]

head(arab.sample)
head(arab.sample.norm)

var.norm <- c()
for (i in 1:n){
  y <- arab.sample.norm[i, ]
  a <- glmmadmb(y ~1, random=  ~1 | group, zeroInflation=F, 
                link="log",  family = "nbinom")
  var.norm[i] <- as.numeric(a$S)  #### output the variance
}

pdf("estimated_variance_normalized_boxplot_hist.pdf")
hist(log(vari), main=paste("number of genes =", n))
hist(log(var.norm), main=paste("number of genes =", n, "normalized"))
boxplot(log(vari), log(var.norm) )
axis(1, 1:2, c("regular","normalized" ), las=1)
dev.off()