
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

### prepare NB data
library(NBPSeq)


norm.factors <- estimate.norm.factors(arab.filter)
lib.size <- colSums(arab.filter)
arab.prepare <- prepare.nb.data(arab.filter, norm.factors=norm.factors)
names(arab.prepare)

#  number of genes to be sampled
n <- 1000
set.seed(120) 


## generate the sequence of genes to be sampled
id <- sample (1:dim(arab.filter)[1], n)  

arab.sample <- arab.filter[id, ]

##  store the variances for gene sampled 
vari <- c()



  for (i in 1:n)
{
  y <- arab.sample[i, ]
  
  a <- glmmadmb(y ~1 + offset(log(lib.size)), random=  ~1 | group,
                zeroInflation=F, link="log",  family = "nbinom")
  
  vari[i] <- as.numeric(a$S)  #### output the variance
}

## considering normalization factors

var.norm <- c();
for (i in 1:n)
{
  y <- arab.sample[i, ]
  
  a <- glmmadmb(y ~1 + offset(log(arab.prepare$eff.lib.sizes)), 
                random=  ~1 | group,
                zeroInflation=F, link="log",  family = "nbinom")
  
  var.norm[i] <- as.numeric(a$S)  #### output the variance
  
}
length(var.norm)
y <- arab.sample[314, ]

a <- glmmadmb(y ~1 + offset(log(arab.prepare$eff.lib.sizes)), 
              random=  ~1 | group,
              zeroInflation=F, link="log",  family = "nbinom")



boxplot(log(vari), log(var.norm))

###  some warning message would pop up. 
#  1: In glmmadmb(y ~ 1, random = ~1 | group, zeroInflation = F, link = "log",  :
# Convergence failed:log-likelihood of gradient= -0.00148958
# convergence issue 

# try y <- c(34, 43, 42, 7, 3, 0, 0, 0, 0, 0, 0 )
# or y <- arab.filter[21251,] 


## draw a histogram of the log variances.
hist(log(vari), main=paste("number of genes =", n))
boxplot(log(vari), log(var.norm))


id2 <- which(log(var.norm)< -10)
## to see which genes have estimated variance of 0
y <- arab.filter[id[id2],]


# we can have glmmPQL estimate the variance, but it cannot be output.
y <- arab.filter[id[25],]

a <- glmmadmb(y~1+offset(log(lib.size)), random=  ~1 | group, zeroInflation=F, 
              link="log",  family = "nbinom")
b <- glmmPQL(y ~1 + offset(off), random=  ~1 | group, family = quasipoisson)

summary(b)
summary(a)
names(summary(b))

# but what is the quantity under (intercetp) of Random effects part
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


pdf("estimated_variance_normalized_boxplot_hist.pdf")
hist(log(vari), main=paste("number of genes =", n))
hist(log(var.norm), main=paste("number of genes =", n, "normalized"))
boxplot(log(vari), log(var.norm) )
axis(1, 1:2, c("regular","normalized" ), las=1)
dev.off()