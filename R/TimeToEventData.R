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

coxme_pval(exp.su,exp,boot=1000)
coxme_icc_ci(exp.su)

###### Mexican Jay problem-solving survival analysis #####
#Latency to solve a door on a puzzle box
jsolv = read.csv("jaySolveData.csv")
jsolv = jsolv[,c(2:6,9)]
colnames(jsolv)[6] = "Time"
jsolv$olre = factor(1:68)

solv.su = coxme(Surv(Time, Solve)~Treatment + (1|ID), data=jsolv)
summary(solv.su)
coxme_pval(solv.su,jsolv,boot = 1000)

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
df.a = read.csv("MEJAintervalAttempts.csv")

Attm.fit = coxme(Surv(Time1, Time2, Attm)~ observe + (1|ID), data=df.a)
summary(Attm.fit)
#random effect variance = 0.0004

# What if the time scale is in minutes instead of seconds?
df.a$newTime1 = df.a$Time1/60
df.a$newTime2 = df.a$Time2/60

Attm.fit2 = coxme(Surv(newTime1, newTime2, Attm)~ observe + (1|ID), data=df.a)
summary(Attm.fit2)
#random effect and fixed effect estimates are the same


## Time to first solve model
df.s = read.csv("timeintervalMEJA.csv")

succ.fit = coxme(Surv(Time1, Time2, success)~observe + (1|ID), data=df.s)
summary(succ.fit)
#random effect variance = 1.277

# If time is minutes instead of seconds
df.s$newTime1 = df.s$Time1/60
df.s$newTime2 = df.s$Time2/60
succ.fit2 = coxme(Surv(newTime1, newTime2, success)~observe + (1|ID), data=df.s)
summary(succ.fit2)
#random effect and fixed effect estimates are the same


## PWE model
df.s$tar = df.s$Time2 - df.s$Time1
interval = data.frame(unique(df.s$Time2))
interval$interval = factor(1:133)
colnames(interval)[1] = "Time2"
jay.int = merge(df.s, interval, by = "Time2", all = T)
meja.pwe.fit = glmer(success ~ observe + interval + (1|ID), data=jay.int,
                family = poisson, offset = log(tar),nAGQ=7,
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#started 6:17pm; ended 4:00pm
summary(meja.pwe.fit)
#random effect variance is essentially zero, but model not converging.


#### Mexican Jay social learning data, ignoring time-varying nature of "Observe" covariate ####
data.naive$deltaT = data.naive$Tot.end - data.naive$Tot.start
Time = aggregate(deltaT ~ ID + Date, data = data.naive, FUN = "sum")
Observe = aggregate(observe ~ ID + Date, data = data.naive, FUN = "sum")
Attempt = aggregate(Attm ~ ID + Date, data = data.naive, FUN = "sum")
Solve = aggregate(success ~ ID + Date, data = data.naive, FUN = "sum")

df_list <- list(Time, Observe, Attempt, Solve) 
SLsumData = df_list %>% reduce(full_join, by=c('ID','Date'))

SLsumData$AttStatus = ifelse(SLsumData$Attm > 0, 1, 0)
SLsumData$SolvStatus = ifelse(SLsumData$success > 0, 1, 0)

## Survival model for attempts with observe as a covariate.
succ.fit2 = coxme(Surv(deltaT, SolvStatus)~observe + (1|ID), data=SLsumData)
summary(succ.fit2)
var(ranef(succ.fit2)$ID)
SLsumData$newT = SLsumData$deltaT/60
succ.fit3 = coxme(Surv(newT, SolvStatus)~observe + (1|ID), data=SLsumData)
summary(succ.fit3) #same results


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
# Pezner et al. 2017 https://academic.oup.com/beheco/article/28/1/154/2453511
ctw = read.csv("CTWemergence.csv")
# "HT" variable indicates hiding time, or the latency to emerge; "Whorls" is a visual indicator of age
# 30 worms received 4 trials per day, across 4 days for a total of 16 trials. 

# 2 individuals have NA values, but it is not explained why. 
# I'll assume these are censored (the worm didn't emerge in the trial time) and give ceiling value of one unit after other highest value
ctw$event = ifelse(is.na(ctw$HT),0,1)
ctw$HT[which(is.na(ctw$HT))]<- 375

# also, a lot of data, so restrict to just one trial per day
ctw2 = ctw[which(ctw$Trial_Total == 1 | ctw$Trial_Total == 2 | ctw$Trial_Total == 3 | ctw$Trial_Total == 4),]

ctw.emerg = coxme(Surv(HT, event)~Whorls + (1|Worm_ID), data=ctw)
summary(ctw.emerg)
#random effect variance = 1.671
var(ranef(ctw.emerg)$Worm_ID) # 1.42
var(random.effects(ctw.emerg)$Worm_ID) #1.42

ctw2$olre = factor(1:120)
ctw.emerg.olre = coxme(Surv(HT, event)~Whorls + (1|Worm_ID) + (1|olre), data=ctw2)
summary(ctw.emerg.olre)
#random effect variance (ID) = 1.672
var(ranef(ctw.emerg.olre)$Worm_ID) # 1.42


## ctw intervals
ctw.int = survSplit(Surv(HT,event) ~ Whorls + Trial_Total + Worm_ID, data = ctw, 
                    cut = ctw$HT, start = "tstart",end = "tstop")

ctw.cox.int = coxme(Surv(tstart,tstop, event)~Whorls + (1|Worm_ID), data=ctw.int)
var <- VarCorr(ctw.cox.int)$Worm_ID[[1]] # random effect variance = 0.81
var/(var + pi^2/6) # 0.33


## PWE model
# for the piecewise exponential model, we need a measure of the time-at-risk for the offset.
ctw.int$tar = ctw.int$tstop - ctw.int$tstart
interval = data.frame(unique(ctw.int$tstop))
interval$interval = factor(1:94)
colnames(interval)[1] = "tstop"
ctw.int = merge(ctw.int, interval, by = "tstop", all = T)

ctw.pwe.fit = glmer(event ~ Whorls + interval + (1|Worm_ID), data=ctw.int,
                family = poisson, offset = log(tar),nAGQ=7)
#started at 3:41pm; ended at 3:43pm with failure to converge warning
summary(ctw.pwe.fit)
#random effect variance = 1.31

## Discrete time glmm
ctw.dtsm.fit = glmer(event ~ Whorls + interval + (1|Worm_ID), data=ctw.int,
                 family = binomial(link="cloglog"), nAGQ=7)
#started 5:24pm; ended 6:07
summary(ctw.dtsm.fit)
#random effect variance = 0.78 ... but warnings it failed to converge


ctw.dtsm.fit2 = glmer(event ~ Whorls + interval + (1|Worm_ID), data=ctw_int,
                     family = binomial(link="cloglog"), nAGQ=7)
#started 3:38pm; ended 
summary(ctw.dtsm.fit2)
#random effect variance = 0.87 ... same as cox model


## quantile ctw intervals - where there are an equal number of events in each of 2 intervals
cutpoints_q <- quantile(ctw$HT, probs = seq(0,1,0.25))
#cutpoints_q <- cutpoints_q[-length(cutpoints_q)] # remove the last value

ctw_int_q <- survSplit(Surv(HT, event) ~ Whorls + Trial_Total + Worm_ID,
                            data = ctw,
                            cut = c(14,20,33), 
                            start = "tstart",
                            end = "tstop")
ctw_int_q$interval = NA
ctw_int_q$interval =  ifelse(ctw_int_q$tstart == 0,1,ctw_int_q$interval)
ctw_int_q$interval =  ifelse(ctw_int_q$tstart == 14,2,ctw_int_q$interval)
ctw_int_q$interval =  ifelse(ctw_int_q$tstart == 20,3,ctw_int_q$interval)
ctw_int_q$interval =  ifelse(ctw_int_q$tstart == 33,4,ctw_int_q$interval)
ctw.dtsm.fit3 = glmer(event ~ Whorls + interval + (1|Worm_ID), data=ctw_int_q,
                      family = binomial(link="cloglog"), nAGQ=7)
VarCorr(ctw.dtsm.fit3)$Worm_ID[[1]]
# 4 intervals random effect variance = 0.82 (2 intervals random effect variance = 0.59)

##### Seed dispersal distance ####
data<-read.csv(file="data_seed_pers.csv")
pm<-subset(data, SPP=="PM") # only one species

pmdistmov<-subset(pm,REMOVE==1) # distance the seed is dispersed, only looking at seeds that were removed from the feeding platform
pmdistmov<-subset(pmdistmov, CONS!=1) # remove rows where the seed was consumed close by the feeding platform
pmdistmov$RECOVERED..Y.N.[which(pmdistmov$DIST..MOVED==15)]<-"Y" # one row seems to indicate the cached seed was found 15m from the feeding platform, but the categorical variable for whether it was removed is an NA 

table(pmdistmov$RECOVERED..Y.N.) # how many cached seeds did they recover (i.e., how much censored data is there?)
# 35% of data are censored

dist = pmdistmov[-which(pmdistmov$ID == "UNK"),c(1,5,6,14,17)] #simplify data frame

# the Recovered column is analogous to the event variable in surivival analysis. Modify it to be an integer
dist$event = ifelse(dist$RECOVERED..Y.N.== "Y",1,0)
quantile(dist$DIST..MOVED, probs = seq(0,1,0.25),na.rm=T)

psych::describe(dist$DIST..MOVED)
dist$DIST..MOVED[which(is.na(dist$DIST..MOVED))]<-1038 # give ceiling value to NAs

dist_int = survSplit(Surv(DIST..MOVED, event) ~ TRT + GRID + ID, data = dist, 
                     cut = c(50,146,331.5), start = "tstart",end = "tstop")
dist.cox.int = coxme(Surv(tstart,tstop, event) ~ 1 + (1|ID), data=dist_int)
VarCorr(dist.cox.int)$ID[[1]] 
# Random effect variance = 0.09


