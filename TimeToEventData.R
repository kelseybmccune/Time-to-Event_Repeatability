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


###### Mexican Jay social learning (includes time-varying covariate) ######
# Performance of naive jays (never interacted with the puzzle box task before) as a function of the number of interactions 
# they observed group members make at the foraging task. There are censored data because some individuals never attempt or solve the task
# 17 naive jays in this sample

SLraw = read.csv("MEJAsocLearnRaw.csv")
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


library(tidyverse)
data.naive = SLraw[which(SLraw$Experienced == "No"),] #remove Experienced jay data
data.naive = data.naive[-which(str_detect(data.naive$Behav, "scrounge")),] #remove scrounges

tmp = data.naive %>% 
  dplyr::select(observe, Attm, success, ID, Group,Trial, Tot.end) %>% 
  arrange(Group, ID, Trial,Tot.end)

### create a dataframe with the time intervals between each event 
#i.e., any time an event happens where a jay did something or saw something, 
#create a time interval and populate rows for each individual for that time interval.
df.tmp = expand.grid('ID' = unique(tmp$ID), 
                     'Time1' = unique(c(tmp$Tot.end, sum.trial.times$Clip.duration.sec.))) %>% 
  arrange(ID, Time1) %>% 
  rename(Time2 = Time1) %>% 
  group_by(ID) %>% 
  mutate(Time1 = lag(Time2)) %>% 
  ungroup() %>% 
  dplyr::select(ID, Time1, Time2) %>% 
  mutate(Time1 = ifelse(is.na(Time1), 0, Time1)) %>% 
  mutate(observe = 0, Attm = 0, success = 0)

groups = sum.data19[,c(1,3)]
df.tmp = merge(df.tmp, groups, by = "ID") #add back in social group ID

for(i in unique(df.tmp$ID)){ #for each bird
  for(j in unique(df.tmp$Time2)){ #for each interval
    tmp2 = data.naive %>% #in original data frame
      filter(ID==i & Tot.end < j & observe == 1) #find the first time the bird observed something
    df.tmp$observe[which(df.tmp$ID==i & df.tmp$Time2==j)] = nrow(tmp2)
    
    tmp2 = data.naive %>% 
      filter(ID==i & Tot.end <= j & Attm == 1) #find the first time the bird attempted
    df.tmp$Attm[which(df.tmp$ID==i & df.tmp$Time2==j)] = nrow(tmp2)
    
    tmp2 = data.naive %>% 
      filter(ID==i & Tot.end <= j & success == 1) #find the first time the bird succeeded
    df.tmp$success[which(df.tmp$ID==i & df.tmp$Time2==j)] = nrow(tmp2)
  }
}

df.tmp$ID = as.factor(as.character(df.tmp$ID))

### Survival analysis with only naive jays and a time-varying covariate of observed interactions by group members
## Attempts model
df.tmp.a = df.tmp[-which(df.tmp$Attm > 1),] 
# only want the data up until the first attempt
df.a = df.tmp.a %>% 
  arrange(ID, Time1) %>% 
  group_by(ID) %>% 
  mutate(testVal = lag(Attm)) %>% 
  ungroup() %>% 
  mutate(testVal = ifelse(is.na(testVal), 0, testVal)) %>% 
  filter(testVal < 1)
# this code insures we are only looking at the time points and number of observations
# up until the jay's first attempt

Attm.fit = coxme(Surv(Time1, Time2, Attm)~ observe + (1|ID), data=df.a)
summary(Attm.fit)
# p = 0.15 - No relationship between observing interactions and the latency to make an attempt
cox.zph(Attm.fit) #proportional hazards assumption not violated

## Calculate 95% confidence intervals for coefficient and hazard ratio
confint(Attm.fit)
exp(confint(Attm.fit))

## Solves model
df.tmp.s = df.tmp[-which(df.tmp$success > 1),] 
df.s = df.tmp.s %>% 
  arrange(ID, Time1) %>% 
  group_by(ID) %>% 
  mutate(testVal = lag(success)) %>% 
  ungroup() %>% 
  mutate(testVal = ifelse(is.na(testVal), 0, testVal)) %>% 
  filter(testVal < 1) 
# this code insures we are only looking at the time points and number of observations
# up until the jay's first success.

succ.fit = coxme(Surv(Time1, Time2, success)~observe + (1|ID), data=df.s)
summary(succ.fit)
# p = 0.03 - a one-unit increase in observations of group members interacting increases 
# the likelihood of naive jays solving by 82% 
# compared to a jay that never observes interactions
cox.zph(succ.fit) #proportional hazards assumption not violated

## Calculate 95% confidence intervals for estimate and hazard ratio
confint(succ.fit)
exp(confint(succ.fit))


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
cox.zph(emerg.fit) # RMR changes over time and so violates proportional hazard assumption

plot.fit = survfit(Surv(Latency.to.emerge, Emerge)~Sex, data = crick)
ggsurvplot(plot.fit, data = crick, fun = 'event', 
                       risk.table = F, pval = F, palette = c("black","#999999"), 
                       ylab = "Proportion", size = 1.5, legend = c(0.7,0.2), 
                       xlab = "Latency to emerge from shelter (sec)", font.tickslab=c(10,"plain","black"), 
                       font.x=c(14,"plain","black"), 
                       font.y=c(14,"plain","black"), font.legend = c(12, "plain","black"), size=0.5, 
                       break.y.by = 0.25, censor=T)




