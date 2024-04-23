####### Grackle exploration survival analysis ######
#Exploration of novel environment

library(coxme)
library(survival)
exp <- read.csv("grackleExpData.csv")
exp.su = coxme(Surv(LatencyFirstLand, event)~ Condition + FlexGroup + (1|BirdID), data=exp)
summary(exp.su)

exp$olre = factor(1:76)
exp.su.olre = coxme(Surv(LatencyFirstLand, event)~ Condition + FlexGroup + (1|BirdID) + (1|olre), data=exp)
summary(exp.su.olre)

exp.po = glmer(round(LatencyFirstLand) ~ Condition + FlexGroup + (1|BirdID), data = exp, family = "poisson")
summary(exp.po)

exp.po.olre = glmer(round(LatencyFirstLand) ~ Condition + FlexGroup + (1|BirdID) + (1|olre), data = exp, family = "poisson")
summary(exp.po.olre)



###### Mexican Jay problem-solving survival analysis #####
#Latency to solve a door on a puzzle box
jsolv = read.csv("jaySolveData.csv")
jsolv$olre = factor(1:68)

solv.su = coxme(Surv(Adjusted, Solve)~Treatment + (1|ID), data=jsolv)
summary(solv.su)

solv.su.olre = coxme(Surv(Adjusted, Solve)~Treatment + (1|ID) + (1|olre), data=jsolv)
summary(solv.su.olre)

solv.po.olre = glmer(round(Adjusted) ~ Treatment + (1|ID) + (1|olre), data = jsolv, family = "poisson")
summary(solv.po.olre)
