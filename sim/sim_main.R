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
#job <- 1250

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


##### modeling  --------------------------------------------------------


###### coxme -------------
# fit a Cox prop hazards model with a frailty term using the coxme package
fit_coxme <- quietly(function(dat) {
  coxme(Surv(survival_time, event) ~ sex + (1 | cluster), data = dat)
})

# save model output
coxme_mod <- fit_coxme(dat)
coxme_mod_res <- coxme_mod$result
coxme_mod_errors <- coxme_mod$messages
coxme_mod_warnings <- coxme_mod$warnings

# extract model estimates if there is no error or warning
if (length(coxme_mod_warnings)==0 && length(coxme_mod_errors)==0) {
  coxme_mod_res <- coxme_mod$result
  beta_coxme <- fixef(coxme_mod_res)
  var_coxme <- VarCorr(coxme_mod_res)$cluster[[1]]
  ICC_coxme <- fr.lognormal(k, s, var_coxme, what = "tau")
  coxme_mod_errors <- NA
  coxme_mod_warnings <- NA
} else { # set to NA and save warning or error message
  coxme_mod_res <- NA
  beta_coxme <- NA
  var_coxme <- NA
  ICC_coxme <- NA
  coxme_mod_errors <- ifelse(length(coxme_mod_errors)==0,
                              NA, paste(coxme_mod_errors, collapse="; "))
  coxme_mod_warnings <- ifelse(length(coxme_mod_warnings)==0,
                                NA, paste(coxme_mod_warnings, collapse="; "))
}




#### coxph1 (normal) -------------
# fit a Cox prop hazards model with a frailty term using the survival package (assuming normal dist)
fit_coxph1 <- quietly(function(dat) {
  coxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gaussian"), dat)
})

# save model output
coxph1_mod <- fit_coxph1(dat)
coxph1_mod_res <- coxph1_mod$result
coxph1_mod_errors <- coxph1_mod$messages
coxph1_mod_warnings <- coxph1_mod$warnings

# extract model estimates if there is no error or warning
if (length(coxph1_mod_warnings)==0 && length(coxph1_mod_errors)==0) {
  coxph1_mod_res <- coxph1_mod$result
  beta_coxph1 <- coxph1_mod_res$coefficients
  var_coxph_normal <- coxph1_mod_res$history$`frailty(cluster, distribution = "gaussian")`$theta
  ICC_coxph1 <- fr.lognormal(1, sigma2.f, var_coxph_normal, what = "tau")
  coxph1_mod_errors <- NA
  coxph1_mod_warnings <- NA
} else { # set to NA and save warning or error message
  coxph1_mod_res <- NA
  beta_coxph1 <- NA
  var_coxph_normal <- NA
  ICC_coxph1 <- NA
  coxph1_mod_errors <- ifelse(length(coxph1_mod_errors)==0,
                              NA, paste(coxph1_mod_errors, collapse="; "))
  coxph1_mod_warnings <- ifelse(length(coxph1_mod_warnings)==0,
                                NA, paste(coxph1_mod_warnings, collapse="; "))
}




#### coxph2 (gamma) -------------
# fit a Cox prop hazards model with a frailty term using the survival package (assuming gamma dist)

fit_coxph2 <- quietly(function(dat) {
  coxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gamma"), dat)
})

# save model output
coxph2_mod <- fit_coxph2(dat)
coxph2_mod_res <- coxph2_mod$result
coxph2_mod_errors <- coxph2_mod$messages
coxph2_mod_warnings <- coxph2_mod$warnings

# extract model estimates if there is no error or warning
if (length(coxph2_mod_warnings)==0 && length(coxph2_mod_errors)==0) {
  coxph2_mod_res <- coxph2_mod$result
  beta_coxph2 <- coxph2_mod_res$coefficients
  var_coxph_gamma <- coxph2_mod_res$history$`frailty(cluster, distribution = "gamma")`$theta # expected to be different
  ICC_coxph2 <- var_coxph_gamma/(var_coxph_gamma + 2)
  coxph2_mod_errors <- NA
  coxph2_mod_warnings <- NA
} else { # set to NA and save warning or error message
  coxph2_mod_res <- NA
  beta_coxph2 <- NA
  var_coxph_gamma <- NA
  ICC_coxph2 <- NA
  coxph2_mod_errors <- ifelse(length(coxph2_mod_errors)==0,
                              NA, paste(coxph2_mod_errors, collapse="; "))
  coxph2_mod_warnings <- ifelse(length(coxph2_mod_warnings)==0,
                                NA, paste(coxph2_mod_warnings, collapse="; "))
}




#### glmm (binomial) -------------
# fit a discrete-time survival model with a frailty term using the glmer function from the lme4 package
fit_glmm <- quietly(function(dat) {
  glmer(event ~ -1 + t_interval + sex + (1|cluster), 
        data=exploded_dat, family=binomial(link="cloglog"))
})

# save model output
glmm_mod <- fit_glmm(dat)
glmm_mod_res <- glmm_mod$result
glmm_mod_errors <- glmm_mod$messages
glmm_mod_warnings <- glmm_mod$warnings

# extract model estimates if there is no error or warning
if (length(glmm_mod_warnings)==0 && length(glmm_mod_errors)==0) {
  glmm_mod_res <- glmm_mod$result
  beta_glmm <- fixef(glmm_mod_res)["sex"]
  var_glmm <- VarCorr(glmm_mod_res)$cluster[[1]]
  ICC_glmm <- var_glmm/(var_glmm + pi^2/6)
  glmm_mod_errors <- NA
  glmm_mod_warnings <- NA
} else { # set to NA and save warning or error message
  glmm_mod_res <- NA
  beta_glmm <- NA
  var_glmm <- NA
  ICC_glmm <- NA
  glmm_mod_errors <- ifelse(length(glmm_mod_errors)==0,
                            NA, paste(glmm_mod_errors, collapse="; "))
  glmm_mod_warnings <- ifelse(length(glmm_mod_warnings)==0,
                              NA, paste(glmm_mod_warnings, collapse="; "))
}




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
                  beta_est = c(beta_coxme,
                               beta_coxph1,
                               beta_coxph2,
                               beta_glmm),
                  ICC = c(ICC_coxme,
                          ICC_coxph1,
                          ICC_coxph2,
                          ICC_glmm),
                  model_warnings = c(coxme_mod_warnings,
                                     coxph1_mod_warnings,
                                     coxph2_mod_warnings,
                                     glmm_mod_warnings),
                  model_errors = c(coxme_mod_errors,
                                   coxph1_mod_errors,
                                   coxph2_mod_errors,
                                   glmm_mod_errors))


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
