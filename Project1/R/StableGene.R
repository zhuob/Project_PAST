
## load the data
source("/Users/Bin/Disk1/Stat/Research/Project2014/Project1/R/prepare.data.R")

##
library(glmmADMB)

mod <- glmmadmb(y~1 + offset(log(off.set)), random = ~(1|group/trt), 
                zeroInflation =F, link="log", family="nbinom")

library(gamlss.mx)


library(lme4)
library(optimx)
##  random-effects Poisson
y <- as.numeric(stable2[4000, ])
id <- 1:dim(stable2)[2]

### mod1 and mod2 has similar results

## trt nested in group
mod1 <- glmer(y ~ 1  + (1|group/trt) + (1|id),
              offset=(log(off.set)), family = poisson, 
              control=glmerControl(optimizer="bobyqa",boundary.tol=1e-7) # advised by Ben Bolker
            # http://stackoverflow.com/questions/21344555/convergence-error-for-development-version-of-lme4
             # control=glmerControl(optimizer="Nelder_Mead")
           #  control=glmerControl(optimizer="optimx",
            #                      optCtrl=list(method="nlminb"))
           )

## negative binomial regression in lme4 package
mod2 <- glmer.nb(y~ 1  + (1|group/trt)+offset(log(off.set))  )

n <- 100
set.seed(13420)
samp <- sample(1:dim(stable2)[1], n)

v <- v2 <- phi <- c()

for ( i in 1: n){
  y <- as.numeric(stable2[samp[i], ])
  print(i)
  mod <- glmmadmb(y~1 + offset(log(off.set)), random = ~(1|group/trt), 
                  zeroInflation =F, link="log", family="nbinom")
 
  v[i] <- as.numeric(mod$S[1]) # the variance of random effect for group
  v2[i] <- as.numeric(mod$S[2]) # the variance of random effect for trt

  phi[i] <- 1/as.numeric(mod$alpha)
}

hist(v2)
var1 <- var2 <- var3 <- c()
for ( i in 1:n){
  y <- as.numeric(stable2[i, ])
  id <- 1:length(y)
  print(i)
  mod1 <- glmer(y ~ 1  + (1|group) + (1 |trt) + (1|id),
                offset=(log(off.set)), family = poisson)
  
  var1[i] <- as.numeric(summary(mod1)$varcor[1]) # variance for individual random effect
  var2[i] <- as.numeric(summary(mod1)$varcor[2]) # variance for trt effect
  var3[i] <- as.numeric(summary(mod1)$varcor[3]) # variance for group effect
}


var1
hist(var3)



# 1707 16137



