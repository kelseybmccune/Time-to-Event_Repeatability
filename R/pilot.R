# pilot

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
#library(survminer)

# data

dat <- read.csv(here("data", "grackleExpData.csv"))

#head(dat)

names(dat)[1] <- "olre"

str(dat)

#dat$olre = factor(1:nrow(dat))

# modeling

coxme0 <-  coxme(Surv(LatencyFirstLand, event) ~ Condition  + (1|BirdID), data=dat)
summary(coxme0)

mu<- (mean(predict(coxme0, type = "risk")) + 0.5*0.7711779)

mu*(exp(0.7711779) - 1)/(mu*(exp(0.7711779) - 1) + 1)
0.7711779/(0.7711779 + trigamma(mu))

model <- lmer(pred ~ 1 + (1|BirdID), data = dat)
summary(model)


coxme1 <- coxme(Surv(LatencyFirstLand, event)~  (1|BirdID) , data=dat)
summary(coxme1)

predict(coxme1, type = "risk")


pois0 <- glmer(round(LatencyFirstLand) ~ Condition  + (1|BirdID) + (1|olre), data = dat, family = "poisson")
summary(pois0)

#pois1 <- glmer(round(LatencyFirstLand) ~ 1  + (1|BirdID) + (1|olre), data = dat, family = "poisson")
#summary(pois1)

#mu <-exp(mean(predict(pois1)) + 0.5*(1.3077 + 0.7483))

#mu*(exp(0.7483) - 1)/(mu*(exp(1.3077 + 0.7483) - 1) + 1)

pois2 <- glmer(event ~ Condition + factor(LatencyFirstLand) + (1|BirdID) , data = dat, family = "poisson")
summary(pois2)


#######

jsolv = read.csv(here("data", "jaySolveData.csv"))
#jsolv$olre = factor(1:68)

solv.su = coxme(Surv(Adjusted, Solve)~Treatment + (1|ID), data=jsolv)
summary(solv.su)

mu <- mean(predict(solv.su , type = "lp") + 0.5* VarCorr(solv.su)[[1]])

VarCorr(solv.su)[[1]]/(VarCorr(solv.su)[[1]] + trigamma(mu))
#########

ctw <- read.csv(here("data","CTWemergence.csv"))
# "HT" variable indicates hiding time, or the latency to emerge; "Whorls" is a visual indicator of age
# 30 worms received 4 trials per day, across 4 days for a total of 16 trials. 

# 2 individuals have NA values, but it is not explained why. I'll assume these are censored (the worm didn't emerge in the trial time)
ctw$event <- ifelse(is.na(ctw$HT),0,1)

ctw.emerg <- coxme(Surv(HT, event)~Whorls + (1|Worm_ID), data=ctw)
summary(ctw.emerg)

predict(ctw.emerg , type = "risk")

mu <- exp(mean(predict(ctw.emerg , type = "lp")) + 0.5*1.030431)

1.030431/(1.030431 + trigamma(mu))

# criket

##### Cricket time-to-emerge data ####
crick = read.csv(here("data","CricketEmergence.csv"))
crick$Cricket = crick$Cricket %>% str_replace(".*-", "")

# data already includes a status column ("Emerge") and time
emerg.fit = coxme(Surv(Latency.to.emerge, Emerge)~ Sex + Mass + RMR + (1|Cricket), data=crick)
summary(emerg.fit)

s2 <- VarCorr(emerg.fit)[[1]]

mu <- exp(mean(predict(emerg.fit , type = "lp")) + 0.5*s2)

s2/(s2 + trigamma(mu))

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

