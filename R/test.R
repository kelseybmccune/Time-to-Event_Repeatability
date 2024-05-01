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
library(parfm)

g <- function(w, k, s, sigma2) {
  -k * w + exp(w) * s + w ^ 2 /  (2 * sigma2)
}

g1 <- function(w, k, s, sigma2) {
  -k + exp(w) * s + w / sigma2
}

g2 <- function(w, k, s, sigma2) {
  exp(w) * s + 1 / sigma2    
} 

Lapl <- Vectorize(function(s, k, sigma2) {
  # Find wTilde = max(g(w)) so that g'(wTilde; k, s, theta) = 0
  WARN <- getOption("warn")
  options(warn = -1)
  wTilde <- optimize(f = g, c(-1e10, 1e10), maximum = FALSE,
                     k = k, s = s, sigma2 = sigma2)$minimum
  options(warn = WARN)
  
  # Approximate the integral via Laplacian method
  res <- (-1) ^ k * 
    exp(-g(w = wTilde, k = k, s = s, sigma2 = sigma2)) /
    sqrt(sigma2 * g2(w = wTilde, k = k, s = s, sigma2 = sigma2))
  return(res)
}, 's')

intTau <- Vectorize(function(x, intTau.sigma2=sigma2) {
  res <- x * 
    Lapl(s = x, k = 0, sigma2 = intTau.sigma2) *
    Lapl(s = x, k = 2, sigma2 = intTau.sigma2)
  return(res)
}, "x")


tauRes <- function(sigma2 = sigma2) {
  4 * integrate(
  f = intTau, lower = 0, upper = Inf, 
  intTau.sigma2 = sigma2)$value - 1
return(tauRes)
}

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

fit0 <- coxph(Surv(time, event) ~ Z1 + Z2 + frailty(cluster, distribution="gaussian"), dat)
summary(fit0)

fit01 <- coxph(Surv(time, event) ~ Z1 + Z2 + frailty(cluster, distribution="gamma"), dat)
summary(fit01)

0.1257205/(0.1257205 + 2)

mod <- parfm(Surv(time, event) ~ Z1 + Z2, data=dat, cluster="cluster", frailty = "gamma")
mod

mod["theta", "ESTIMATE"]
#tau(mod)
mod1 <- parfm(Surv(time, event) ~ Z1 + Z2, data=dat, cluster="cluster", frailty = "lognormal")
mod1
#fit0e <- coxph(Surv(time*1000, event) ~ Z1 + Z2 + frailty(cluster, distribution="gaussian"), dat)
#summary(fit0e)



predict(fit0, type = "lp")

bhest <- basehaz(fit0)
haz <- exp(diff(bhest[, 1])*diff(bhest[, 2]))


fit0b <- coxph(Surv(time, event) ~ Z1 + Z2 + frailty(cluster, distribution="gamma"), dat, 
              method="breslow", x=TRUE, y=TRUE)
summary(fit0b)

fit1 <- coxme(Surv(time, event) ~ Z1 + Z2 + (1|cluster), dat)
summary(fit1)

fr.lognormal(k,s,0.2384844,what = "tau")
#var(ranef(fit1)$cluster)


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
                  data=ppd, family=poisson)
summary(fit2)
var(pull(ranef(fit2)$cluster))

fit3 <- glmer(m~-1+t+z1+z2+(1|cluster)+offset(log(N)), 
              data=ppd, family=binomial(link="cloglog"))

summary(fit3)
var(pull(ranef(fit3)$cluster))


var <- VarCorr(fit2)$cluster[1,1]
mu <- exp(0.5*var)
mu2 <- mean(ppd$m)

var/(var + trigamma(mu))
var/(var + log(1/mu + 1))
mu*(exp(var) - 1)/(mu*(exp(var) - 1) + 1)


fr.lognormal(k,s,0.2384844,what = "tau")
#var/(var + pi^2/6)
#var/(var + trigamma(mu2))
var/(var + log(1/mu2 + 1))
mu2*(exp(var) - 1)/(mu2*(exp(var) - 1) + 1)
