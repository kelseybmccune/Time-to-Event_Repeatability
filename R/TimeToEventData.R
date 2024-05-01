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
jsolv = jsolv[,c(2:6,9)]
colnames(jsolv)[6] = "Time"
jsolv$olre = factor(1:68)

solv.su = coxme(Surv(Time, Solve)~Treatment + (1|ID), data=jsolv)
summary(solv.su)
var(ranef(solv.su)$ID)

solv.su.olre = coxme(Surv(Time, Solve)~Treatment + (1|ID) + (1|olre), data=jsolv)
summary(solv.su.olre)

solv.po.olre = glmer(round(Time) ~ Treatment + (1|ID) + (1|olre), data = jsolv, family = "poisson")
summary(solv.po.olre)


###### Mexican Jay social learning (includes time-varying covariate) ######
# Performance of naive jays (never interacted with the puzzle box task before) as a function of the number of interactions 
# they observed group members make at the foraging task. There are censored data because some individuals never attempt or solve the task
# 17 naive jays in this sample

# SLraw = read.csv("MEJAsocLearnRaw.csv")
# These data contain a row for every time the focal individual attempted or solved at the novel task, OR anytime they observed a 
# group member interact at the task.

# We want to know the effect of observing group members on naive jays performance 
# Observations is a time-varying covariate because it increases with time
# So we create a dataframe that has a unique row for each individual 
# at each unique Time that something happened
# (i.e., every time someone sees something or does something).
# We will have more intervals than we need for each individual.  
# In other words, two subsequent intervals could have 
# identical traits, but that shouldn't matter.  

## Time to first attempt model first
write.csv(df.a, "MEJAintervalAttempts.csv")
df.a = read.csv("MEJAintervalAttempts.csv")

Attm.fit = coxme(Surv(Time1, Time2, Attm)~ observe + (1|ID), data=df.a)
summary(Attm.fit)
# p = 0.15 - No relationship between observing interactions and the latency to make an attempt
cox.zph(Attm.fit) #proportional hazards assumption not violated

## Calculate 95% confidence intervals for coefficient and hazard ratio
confint(Attm.fit)
exp(confint(Attm.fit))


# What if the time scale is in minutes instead of seconds?
df.a$newTime1 = df.a$Time1/60
df.a$newTime2 = df.a$Time2/60

Attm.fit2 = coxme(Surv(newTime1, newTime2, Attm)~ observe + (1|ID), data=df.a)
summary(Attm.fit2)


## Time to first solve model
df.s = read.csv("timeintervalMEJA.csv")

succ.fit = coxme(Surv(Time1, Time2, success)~observe + (1|ID), data=df.s)
summary(succ.fit)
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




#### Mexican Jay social learning data, ignoring time-varying nature of "Observe" covariate ####
data.naive$deltaT = data.naive$Tot.end - data.naive$Tot.start
Time = aggregate(deltaT ~ ID + Trial, data = data.naive, FUN = "sum")
Observe = aggregate(observe ~ ID + Trial, data = data.naive, FUN = "sum")
Attempt = aggregate(Attm ~ ID + Trial, data = data.naive, FUN = "sum")
Solve = aggregate(success ~ ID + Trial, data = data.naive, FUN = "sum")

df_list <- list(Time, Observe, Attempt, Solve) 
SLsumData = df_list %>% reduce(full_join, by=c('ID','Trial'))

SLsumData$AttStatus = ifelse(SLsumData$Attm > 0, 1, 0)
SLsumData$SolvStatus = ifelse(SLsumData$success > 0, 1, 0)

## Survival model for attempts with observe as a covariate.
succ.fit2 = coxme(Surv(deltaT, SolvStatus)~observe + (1|ID), data=SLsumData)
summary(succ.fit2)

attm.fit2 = coxme(Surv(deltaT, AttStatus)~observe + (1|ID), data=SLsumData)
summary(attm.fit2)


##### Cricket time-to-emerge data ####
crick = read.csv("CricketEmergence.csv")
crick$Cricket = crick$Cricket %>% str_replace(".*-", "")

# data already includes a status column ("Emerge") and time
emerg.fit = coxme(Surv(Latency.to.emerge, Emerge)~Sex + Mass + RMR + (1|Cricket), data=crick)
summary(emerg.fit)
var(ranef(emerg.fit)$Cricket)
cox.zph(emerg.fit) # RMR changes over time and so violates proportional hazard assumption

cph = coxph(Surv(Latency.to.emerge, Emerge)~Sex + Mass + RMR + frailty(Cricket, distribution = "gaussian"), data=crick)
summary(cph)


plot.fit = survfit(Surv(Latency.to.emerge, Emerge)~Sex, data = crick)
ggsurvplot(plot.fit, data = crick, fun = 'event', 
                       risk.table = F, pval = F, palette = c("black","#999999"), 
                       ylab = "Proportion", size = 1.5, legend = c(0.7,0.2), 
                       xlab = "Latency to emerge from shelter (sec)", font.tickslab=c(10,"plain","black"), 
                       font.x=c(14,"plain","black"), 
                       font.y=c(14,"plain","black"), font.legend = c(12, "plain","black"), size=0.5, 
                       break.y.by = 0.25, censor=T)

# define intervals in data using survSplit function for PWE and Discrete time survival models
# I chose unequal interval sizes so that there was one event per interval
crick.int = survSplit(Surv(Latency.to.emerge,Emerge) ~ Sex + Mass + RMR + Cricket, data = crick, 
                      cut = crick$Latency.to.emerge, start = "tstart",end = "tstop")

emerg.int.fit = coxme(Surv(tstart, tstop, Emerge)~Sex + Mass + RMR + (1|Cricket), data=crick.int)
summary(emerg.int.fit) #all estimates remain the same as in the non-interval cox model.
# random effect variance = 8.210121e-05

# for the piecewise exponential model, we need a measure of the time-at-risk for the offset.
crick.int$tar = crick.int$tstop - crick.int$tstart
interval = data.frame(unique(crick.int$tstop))
interval$interval = factor(1:121)
colnames(interval)[1] = "tstop"
crick.int = merge(crick.int, interval, by = "tstop", all = T)

pwe.fit = glmer(Emerge ~ Sex + Mass + RMR + interval + (1|Cricket), data=crick.int,
                family = poisson, offset = log(tar),nAGQ=7)
#started at 10:57am; ended at 11:58 with failure to converge warning
summary(pwe.fit)
# random effect variance = 5.423e-05

dtsm.fit = glmer(Emerge ~ Sex + Mass + RMR+ interval + (1|Cricket), data=crick.int,
                 family = binomial(link="cloglog"), nAGQ=7)
#started 12:00pm; ended 1:05pm with failure to converge warning
summary(dtsm.fit)
#random effect variance = 3.614e-05




#### Christmas tree worm time-to-emerge data ####
ctw = read.csv("CTWemergence.csv")
# "HT" variable indicates hiding time, or the latency to emerge; "Whorls" is a visual indicator of age
# 30 worms received 4 trials per day, across 4 days for a total of 16 trials. 

# 2 individuals have NA values, but it is not explained why. 
# I'll assume these are censored (the worm didn't emerge in the trial time) and give ceiling value of one unit after other highest value
ctw$event = ifelse(is.na(ctw$HT),0,1)
ctw$HT[which(is.na(ctw$HT))]<- 375

# also, a lot of data, so restrict to just one trial per day
ctw2 = ctw[which(ctw$Trial_Total == 1 | ctw$Trial_Total == 2 | ctw$Trial_Total == 3 | ctw$Trial_Total == 4),]

ctw.emerg = coxme(Surv(HT, event)~Whorls + (1|Worm_ID), data=ctw2)
summary(ctw.emerg)
#random effect variance = 1.671
var(ranef(ctw.emerg)$Worm_ID) # 1.42

ctw2$olre = factor(1:120)
ctw.emerg.olre = coxme(Surv(HT, event)~Whorls + (1|Worm_ID) + (1|olre), data=ctw2)
summary(ctw.emerg.olre)
#random effect variance (ID) = 1.672
var(ranef(ctw.emerg.olre)$Worm_ID) # 1.42


## ctw intervals
ctw.int = survSplit(Surv(HT,event) ~ Whorls + Trial_Total + Worm_ID, data = ctw2, 
                    cut = ctw2$HT, start = "tstart",end = "tstop")


ctw.cox.int = coxme(Surv(tstart,tstop, event)~Whorls + (1|Worm_ID), data=ctw.int)
summary(ctw.cox.int)
#random effect variance = 1.671 ... same as cox model with continuous data
var(ranef(ctw.cox.int)$Worm_ID) # 1.42


## PWE model
# for the piecewise exponential model, we need a measure of the time-at-risk for the offset.
ctw.int$tar = ctw.int$tstop - ctw.int$tstart
interval = data.frame(unique(ctw.int$tstop))
interval$interval = factor(1:46)
colnames(interval)[1] = "tstop"
ctw.int = merge(ctw.int, interval, by = "tstop", all = T)

ctw.pwe.fit = glmer(event ~ Whorls + interval + (1|Worm_ID), data=ctw.int,
                family = poisson, offset = log(tar),nAGQ=7)
#started at 3:41pm; ended at 3:43pm with failure to converge warning
summary(ctw.pwe.fit)
#random effect variance = 1.31


ctw.dtsm.fit = glmer(event ~ Whorls + interval + (1|Worm_ID), data=ctw.int,
                 family = binomial(link="cloglog"), nAGQ=7)
#started 3:38pm; ended 
summary(ctw.dtsm.fit)
#random effect variance = 1.67 ... same as cox model


