#################################################################
# Simulation studies of time-to-event models with repeatability #
#################################################################

# load libraries
library(pacman) # checks if package is installed, if not installs it
p_load(survival, coxme, lme4, here, dplyr, purrr)

# load function
source("R/function.R")

##### load sim conditions ----------------------------------------------

# load job array table
tab <- read.csv("sim/job_array.csv")

# get job number from pbs script
job <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

# get current job information
name <- tab$name[tab$job_number == job] 
scen <- tab$scenario[tab$job_number == job] 
seed <- tab$sim[tab$job_number == job]

# load parameter conditions for current job
n <- tab$n[job]                           # the number of observations
n_clusters <- tab$n_clusters[job]         # the number of individuals (cluster)
sigma2.f <- tab$sigma2.f[job]             # variance of frailty random effect
censoring_prop <- tab$censoring_prop[job] # proportion of censored data
intervals <- tab$intervals[job]           # number of time intervals

# set seed for reproducibility
set.seed(seed)

# number of repetitions per individual (cluster)
repetitions <- n/n_clusters


##### data generation  ----------------------------------------------

# generate covariate (sex)
sex <- rbinom(n, size = 1, prob = 0.5)[rep(1:n_clusters, each = repetitions)]
# set true parameter value
beta <- 2

# generate cluster id variable
cluster <- rep(1:n_clusters, each = repetitions)

# generate frailty term (random effect) at the cluster level
frailty <- rnorm(n_clusters, sd = sqrt(2))[cluster]

# generate survival times from an exponential distribution
hazard <- exp(beta * sex + frailty)
survival_time <- rexp(n, hazard)

# obtain censoring time 
censoring_time <- quantile(survival_time, 1-censoring_prop)
event <- ifelse(survival_time <= censoring_time, 1, 0)
# set survival time to censoring time if it is less than the censoring time
survival_time <- pmin(survival_time, censoring_time) 

#### create a data frame (dat)
dat <- data.frame(observation = 1:n,
                  cluster = cluster,
                  sex = sex,
                  survival_time = survival_time,
                  event = event)


#### create a data frame (exploded_dat)
exploded_dat <- survSplit(Surv(survival_time, event) ~ sex + cluster, data = dat,
                          cut = (1:intervals)*(max(dat$survival_time)/intervals), 
                          start = "tstart",end = "tstop",  zero=0, id = "observation")
# add time interval
exploded_dat$t_interval <- as.factor(exploded_dat$tstart)


##### modeling  ----------------------------------------------

##### TODO
### - wrap models into a function and use quietly
### - store warning message if there is one
### - use if statement to extract model estimates or set to zero if there is a warning message



# Fit a Cox proportional hazards model with a frailty term using the coxme package
coxme_fit <- function(data=dat){
  coxme(Surv(survival_time, event) ~ sex + (1|cluster), data = dat)
}


# Fit a Cox proportional hazards model with a frailty term using the survival package (assuming normal dist)
coxph_fit1 <- try(coxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gaussian"), dat))


# Fit a Cox proportional hazards model with a frailty term using the survival package (assuming gamma dist)
coxph_fit2 <- try(coxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gamma"), dat))


# Fit a discrete-time survival model with a frailty term using the glmer function from the lme4 package
glmm_fit <- try(glmer(event~ -1 + t_interval + sex + (1|cluster), 
                  data=exploded_dat, family=binomial(link="cloglog")))



# # set up quietly to extract error and warning messages
qcoxme <- quietly(coxme)
qcoxph <- quietly(coxph)
qglmer <- quietly(glmer)
qcoxme_fit <- qcoxme(Surv(survival_time, event) ~ sex + (1|cluster), data = dat)
qcoxph_fit1 <- qcoxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gaussian"), dat)
qcoxph_fit2 <- qcoxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gamma"), dat)
qglmm_fit <- qglmer(event~ -1 + t_interval + sex + (1|cluster), data=exploded_dat, family=binomial(link="cloglog"))


##### extract model estimates  ---------------------------------

var_coxme <- VarCorr(coxme_fit)$cluster[[1]]
var_coxph_normal <- coxph_fit1$history$`frailty(cluster, distribution = "gaussian")`$theta
var_coxph_gamma <- coxph_fit2$history$`frailty(cluster, distribution = "gamma")`$theta # expected to be different
var_glmm <- VarCorr(glmm_fit)$cluster[[1]]


##### save results  ----------------------------------------------

res <- data.frame(name = rep(name, 4),
                  scenario = rep(scen, 4),
                  simseed_id = rep(seed, 4),
                  n = rep(n, 4),
                  n_clusters = rep(n_clusters, 4),
                  sigma2.f = rep(sigma2.f, 4),
                  censoring_prop = rep(censoring_prop, 4),
                  intervals = rep(intervals, 4),
                  model = c("coxme",
                            "coxph_normal",
                            "coxph_gamma",
                            "glmm"),
                  sigma2.f_est = c(var_coxme,
                                   var_coxph_normal,
                                   var_coxph_gamma,
                                   var_glmm),
                  beta_est = c(fixef(coxme_fit),
                               coxph_fit1$coefficients,
                               coxph_fit2$coefficients,
                               fixef(glmm_fit)["sex"]),
                  ICC = c(fr.lognormal(1, sigma2.f, var_coxme, what = "tau"),
                          fr.lognormal(1, sigma2.f, var_coxph_normal, what = "tau"),
                          var_coxph_gamma/(var_coxph_gamma + 2),
                          var_glmm/(var_glmm + pi^2/6)))


# save the output according to the job array:
save(list = "res", file = paste0("results/raw/res_", job, ".RDATA"))

# # fixed effects
# fixef(coxme_fit)
# fixef(glmm_fit)["sex"]
# coxph_fit1$coefficients
# coxph_fit2$coefficients
# 
# # ICC tau and glmmer
# fr.lognormal(k,s,var_coxme, what = "tau")
# fr.lognormal(k,s,var_coxph_normal, what = "tau")
# var_glmm/(var_glmm + pi^2/6) # expected to be different
# var_coxph_gamma/(var_coxph_gamma + 2)
