# test

remotes::install_github("mcdonohue/phmm")


library(phmm)
library(coxme)
library(here)
library(survival)
library(lme4)
library(lmerTest)
library(tidyverse)
library(Matrix)
library(eha)
## Loading required package:  survival
## Loading required package:  lattice
## Loading required package:  Matrix
n <- 50 # total sample size
nclust <- 5 # number of clusters 
clusters <- rep(1:nclust,each=n/nclust) 
beta0 <- c(1,2)
set.seed(13)
Z <-cbind(Z1=sample(0:1,n,replace=TRUE), Z2=sample(0:1,n,replace=TRUE), Z3=sample(0:1,n,replace=TRUE))
b <- cbind(rep(rnorm(nclust), each=n/nclust), rep(rnorm(nclust), each=n/nclust))
Wb <- matrix(0,n,2)
for( j in 1:2) Wb[,j] <- Z[,j]*b[,j]
Wb <- apply(Wb,1,sum)
T <- -log(runif(n,0,1))*exp(-Z[,c('Z1','Z2')]%*%beta0-Wb) 
C <- runif(n,0,1)
time <- ifelse(T<C,T,C)
event <- ifelse(T <= C,1,0)
sum(event)

dat <- data.frame(Z) 
dat$cluster <- clusters 
dat$time <- time 
dat$event <- event

fit0 <- coxph(Surv(time, event) ~ Z1 + Z2 + frailty(cluster), dat, 
              method="breslow", x=TRUE, y=TRUE)
summary(fit0)

fit1 <- coxme(Surv(time, event) ~ Z1 + Z2 + (1|cluster), dat)
summary(fit1)
var(ranef(fit1)$cluster)


var <- VarCorr(fit1)$cluster[[1]]
mu <- exp(0.5*var)

var/(var + trigamma(mu))

mu*(exp(var) - 1)/(mu*(exp(var) - 1) + 1)
#VarCorr(fit1)$cluster[[1]]/(VarCorr(fit1)$cluster[[1]] + trigamma(1))

# ppd <- as.data.frame(as.matrix(pseudoPoisPHMM(fit0)))
# # pois likelihood
# poisl <- c()
# eventtimes <- sort(dat$time[dat$event == 1])
# 
# for(h in 1:length(eventtimes)){
#   js <- ppd$time == eventtimes[h] & ppd$m >= 1 # j star j <- ppd$time == eventtimes[h]
#   if(sum(js) > 1) 
#     stop("tied event times")
#   poisl <- c(poisl,
#              ppd[js, "N"]*exp(-1)*exp(ppd[js, "linear.predictors"])/ 
#                sum(ppd[j, "N"]*exp(ppd[j, "linear.predictors"])))
# }
# 
# sum(log(poisl))
# 
# ppd$t <- as.factor(ppd$time)

# phmm

fit.phmm <- phmm(Surv(time, event) ~ Z1 + Z2 + (1|cluster), dat, Gbs = 100, Gbsvar = 1000, VARSTART = 1,
                 NINIT = 10, MAXSTEP = 100, CONVERG=90)
summary(fit.phmm)

# creating data

ppd <- as.data.frame(as.matrix(pseudoPoisPHMM(fit.phmm)))


# Poisson GLMM

ppd$t <- as.factor(ppd$time) 
fit2 <- glmer(m~-1+t+z1+z2+(1|cluster)+offset(log(N)), 
                  data=ppd, family=poisson, nAGQ=0)

summary(fit2)
var(pull(ranef(fit2)$cluster))


var <- VarCorr(fit2)$cluster[1,1]
mu <- exp(mean(predict(fit2)) + 0.5*VarCorr(fit2)$cluster[1,1])

var/(var + trigamma(mu))
mu*(exp(var) - 1)/(mu*(exp(var) - 1) + 1)

