# pilot

# check - the veriance is the same in cloglog models..... 

# some thoughts

# Base rate can be obtained by Poisson expansion and get a rate of event aross all the data
# Variance from comxe can be used to get S2b 


#options(repos = c(CRAN = "https://cloud.r-project.org"))
#utils::install.packages("Matrix")
#utils::install.packages("lme4")

# library

library(coxme)
library(here)
library(survival)
library(lme4)
library(lmerTest)
library(tidyverse)
library(Matrix)
library(eha)
library(phmm)
#library(survminer)

# data

dat <- read.csv(here("data", "grackleExpData.csv"))

#head(dat)

names(dat)[1] <- "olre"

str(dat)

#dat$olre = factor(1:nrow(dat))

# modeling

coxme0 <-  coxme(Surv(LatencyFirstLand, event) ~ 1  + (1|BirdID), data=dat)
#summary(coxme0)
var <- VarCorr(coxme0)$BirdID[[1]]
fr.lognormal(k,s,var,what = "tau")
var/(var + pi^2/6)

mu <- mean(predict(coxme0, type = "risk"))

var/(var + trigamma(mu))
var/(var + log(1/mu + 1))
var/(var + 1/mu)
mu*(exp(var) - 1)/(mu*(exp(var) - 1) + 1)

dat$BirdID <- as.numeric(as.factor(dat$BirdID))
mod <- parfm(Surv(LatencyFirstLand, event) ~ Condition, data=dat, cluster="BirdID", frailty = "gamma", method = "BFGS")
mod
0.532/(0.532 + 2)

dat$BirdID <- as.numeric(as.factor(dat$BirdID))
fit.phmm <- phmm(Surv(LatencyFirstLand, event) ~ Condition + (1|BirdID), dat,Gbs = 100, Gbsvar = 1000, VARSTART = 1,
                 NINIT = 10, MAXSTEP = 100, CONVERG=90)
summary(fit.phmm)

ppd <- as.data.frame(as.matrix(pseudoPoisPHMM(fit.phmm)))

mu2 <- mean(ppd$m)

#basehaz(coxme0)

var <- VarCorr(coxme0)$BirdID[[1]]

mu<- (0.5*var)

mu*(exp(var) - 1)/(mu*(exp(var) - 1) + 1)
var/(var + trigamma(mu))
var/(var + log(1/mu + 1))
var/(var + 1/mu)

mu2*(exp(var) - 1)/(mu2*(exp(var) - 1) + 1)
var/(var + trigamma(mu2))
var/(var + log(1/mu2 + 1))
var/(var + 1/mu2)

#model <- lmer(pred ~ 1 + (1|BirdID), data = dat)
#summary(model)


coxme1 <- coxme(Surv(LatencyFirstLand, event)~  (1|BirdID) , data=dat)
summary(coxme1)

predict(coxme1, type = "risk")


#pois0 <- glmer(round(LatencyFirstLand) ~ Condition  + (1|BirdID) + (1|olre), data = dat, family = "poisson")
#summary(pois0)

#pois1 <- glmer(round(LatencyFirstLand) ~ 1  + (1|BirdID) + (1|olre), data = dat, family = "poisson")
#summary(pois1)

#mu <-exp(mean(predict(pois1)) + 0.5*(1.3077 + 0.7483))

#mu*(exp(0.7483) - 1)/(mu*(exp(1.3077 + 0.7483) - 1) + 1)

pois2 <- glmer(event ~ Condition + factor(LatencyFirstLand) + (1|BirdID) , data = dat, family = "poisson")
summary(pois2)


#######

jsolv = read.csv(here("data", "jaySolveData.csv"))
#jsolv$olre = factor(1:68)

solv.su = coxme(Surv(Adjusted, Solve)~1 + (1|ID), data=jsolv)
#summary(solv.su)


jsolv$ID <- as.numeric(as.factor(jsolv$ID))
#mod <- parfm(Surv(Adjusted, Solve) ~ Treatment, data=jsolv, cluster="ID", frailty = "lognormal", method="BFGS")
#mod

var <- VarCorr(solv.su)$ID[[1]]
fr.lognormal(k,s,var,what = "tau")
var/(var + pi^2/6)

mu <- mean(predict(solv.su, type = "risk"))

var/(var + trigamma(mu))
var/(var + log(1/mu + 1))
var/(var + 1/mu)


##########
jsolv$ID <- as.numeric(as.factor(jsolv$ID))
fit.phmm <- phmm(Surv(Adjusted, Solve) ~ Treatment + (1|ID), jsolv,
                 Gbs = 100, Gbsvar = 1000, VARSTART = 1,
                 NINIT = 10, MAXSTEP = 100, CONVERG=90)
summary(fit.phmm)

ppd <- as.data.frame(as.matrix(pseudoPoisPHMM(fit.phmm)))

ppd$t <- as.factor(ppd$time) 
fit2 <- glmer(m~-1+t+z1+(1|cluster)+offset(log(N)), 
              data=ppd, family=poisson)
summary(fit2)


mean(ppd$m)

var <- VarCorr(solv.su)$ID[[1]]

mu <- mean(0.5* var)
mu2 <- mean(ppd$m)

# mu*(exp(var) - 1)/(mu*(exp(var) - 1) + 1)
# var/(var + trigamma(mu))
# var/(var + log(1/mu + 1))
# var/(var + 1/mu)
fr.lognormal(k,s,var,what = "tau")
mu2*(exp(var) - 1)/(mu2*(exp(var) - 1) + 1)
#var/(var + trigamma(mu2))
var/(var + log(1/mu2 + 1))
#var/(var + 1/mu2)

#########

ctw <- read.csv(here("data","CTWemergence.csv"))
# "HT" variable indicates hiding time, or the latency to emerge; "Whorls" is a visual indicator of age
# 30 worms received 4 trials per day, across 4 days for a total of 16 trials. 

# 2 individuals have NA values, but it is not explained why. I'll assume these are censored (the worm didn't emerge in the trial time)
ctw$event = ifelse(is.na(ctw$HT),0,1)
ctw$HT[which(is.na(ctw$HT))]<- 375

# also, a lot of data, so restrict to just one trial per day
ctw2 = ctw[which(ctw$Trial_Total == 1 | ctw$Trial_Total == 2 | ctw$Trial_Total == 3 | ctw$Trial_Total == 4),]

ctw.emerg <- coxme(Surv(HT, event)~1 + (1|Worm_ID), data=ctw)
#summary(ctw.emerg)
var <- VarCorr(ctw.emerg)$Worm_ID[[1]]
fr.lognormal(k,s,var,what = "tau")
var/(var + pi^2/6)

mu <- mean(predict(ctw.emerg, type = "risk"))

var/(var + trigamma(mu))
var/(var + log(1/mu + 1))
var/(var + 1/mu)
mu*(exp(var) - 1)/(mu*(exp(var) - 1) + 1)

ctw$Worm_ID <- as.numeric(as.factor(ctw$Worm_ID))
#mod <- parfm(Surv(HT, event) ~ Whorls, data=ctw, cluster="Worm_ID", frailty = "gamma")
#mod


ctw$Worm_ID <- as.numeric(as.factor(ctw$Worm_ID))
fit.phmm <- phmm(Surv(HT, event) ~ Whorls + (1|Worm_ID), ctw,
                 Gbs = 100, Gbsvar = 1000, VARSTART = 1,
                 NINIT = 10, MAXSTEP = 100, CONVERG=90)
summary(fit.phmm)

mu2 <- mean(ppd$m)

predict(ctw.emerg , type = "risk")

var <- VarCorr(ctw.emerg)$Worm_ID[[1]]

mu <- mean(predict(coxme0, type = "risk"))

var/(var + trigamma(mu))
var/(var + log(1/mu + 1))
var/(var + 1/mu)
mu*(exp(var) - 1)/(mu*(exp(var) - 1) + 1)

mu2*(exp(var) - 1)/(mu2*(exp(var) - 1) + 1)
var/(var + trigamma(mu2))
var/(var + log(1/mu2 + 1))
fr.lognormal(k,s,var,what = "tau")
var/(var + 1/mu2)

# exploding

# also, a lot of data, so restrict to just one trial per day
#ctw2 = ctw[which(ctw$Trial_Total == 1 | ctw$Trial_Total == 2 | ctw$Trial_Total == 3 | ctw$Trial_Total == 4),]

ctw.emerg <- coxme(Surv(HT, event)~Whorls + (1|Worm_ID), data=ctw2)
summary(ctw.emerg)
var <- VarCorr(ctw.emerg)$Worm_ID[[1]]
fr.lognormal(k,s,var,what = "tau")

ctw.emerg2 <- coxph(Surv(HT, event)~Whorls + frailty(Worm_ID, distribution="gamma"), data=ctw2)
summary(ctw.emerg2)

var2 <- ctw.emerg2$history$`frailty(Worm_ID, distribution = "gamma")`$theta
var2/(var2 + 2)

## ctw intervals
ctw.int = survSplit(Surv(HT,event) ~ Whorls + Trial_Total + Worm_ID, data = ctw2, 
                    cut = ctw2$HT, start = "tstart",end = "tstop")

mu2 <- mean(ctw.int$event)

#mu <- exp(mean(predict(emerg.fit , type = "lp")) + 0.5*var)

mu2*(exp(var) - 1)/(mu2*(exp(var) - 1) + 1)
var/(var + log(1/mu2 + 1))



# criket

##### Cricket time-to-emerge data ####
crick = read.csv(here("data","CricketEmergence.csv"))
crick$Cricket = crick$Cricket %>% str_replace(".*-", "")

# data already includes a status column ("Emerge") and time
emerg.fit = coxme(Surv(Latency.to.emerge, Emerge)~ Sex + Mass + RMR + (1|Cricket), data=crick)
summary(emerg.fit)

var <- as.numeric(VarCorr(emerg.fit)[[1]])
fr.lognormal(k,s,var,what = "tau")
var/(var + pi^2/6)


mu <- exp(mean(predict(emerg.fit , type = "lp")) + 0.5*var)

mu*(exp(var) - 1)/(mu*(exp(var) - 1) + 1)
var/(var + trigamma(mu))
var/(var + log(1/mu + 1))
var/(var + 1/mu)

#####

## Time to first attempt model first


## Time to first solve model
df.s <- read.csv(here("data", "timeintervalMEJA.csv"))

succ.fit <- coxme(Surv(Time1, Time2, success)~observe + (1|ID), data=df.s)
summary(succ.fit)

fit0 <- coxph(Surv(Time1, Time2, success) ~ observe + frailty(ID, distribution="gaussian"), df.s)
summary(fit0)

var <- 1.179973
fr.lognormal(k,s,var ,what = "tau")
var/(var + pi^2/6)

fit01 <- coxph(Surv(Time1, Time2, success) ~ observe+ frailty(ID, distribution="gamma"), df.s)
summary(fit01)

#5e-07/(5e-07 + 2)
df.s$ID <- as.numeric(as.factor(df.s$ID))
mod <- parfm(Surv(Time1, Time2, success) ~ observe, data=df.s, cluster="ID", frailty = "gamma")
mod

mod["theta", "ESTIMATE"]
#tau(mod)
mod1 <- parfm(Surv(Time1, Time2, success) ~ observe, data=df.s, cluster="ID", frailty = "lognormal")
mod1

# does not run!!
df.s$ID <- as.numeric(as.factor(df.s$ID))
fit.phmm <- phmm(Surv(Time1, success) ~ observe + (1|ID), df.s,Gbs = 100, Gbsvar = 1000, VARSTART = 1,
                 NINIT = 10, MAXSTEP = 100, CONVERG=90)
summary(fit.phmm)


ppd <- as.data.frame(as.matrix(pseudoPoisPHMM(fit.phmm)))

ppd$t <- as.factor(ppd$time) 
fit2 <- glmer(m~-1+t+z1+(1|cluster)+offset(log(N)), 
              data=ppd, family=poisson)
summary(fit2)

# p = 0.03 - a one-unit increase in observations of group members interacting increases 
# the likelihood of naive jays solving by 82% 
# compared to a jay that never observes interactions
cox.zph(succ.fit) #proportional hazards assumption not violated

## Calculate 95% confidence intervals for estimate and hazard ratio
confint(succ.fit)
exp(confint(succ.fit))

# If time is minutes instead of seconds
df.s$newTime1 = df.s$Time1/60
df.s$newTime2 = df.s$Time2/60
succ.fit2 = coxme(Surv(newTime1, newTime2, success)~observe + (1|ID), data=df.s)
summary(succ.fit2)




########
#more tests 

dat <- data.frame(enter = rep(0, 4), exit = 1:4,
                  event = rep(1, 4), x = c(0, 1, 0, 1))
dat

fit <- coxreg(Surv(enter, exit, event)~ x, data = dat)
summary(fit)

predict(fit)


datB <- toBinary(dat)
datB

fit2 <- glm(event~ riskset + x, family = poisson,
            data = datB)

summary(fit2)


fit.ml <- coxreg(Surv(enter, exit, event)~ x,
                 method = 'ml', data = dat)

fit.b <- glm(event~ riskset + x,
             family = binomial(link = 'cloglog'),
             data = datB)

summary(fit.ml)
summary(fit.b)

