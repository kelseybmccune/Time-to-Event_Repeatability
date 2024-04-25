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

predict(coxme0, type = "risk")

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

predict(solv.su , type = "risk")





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

