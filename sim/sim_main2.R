#################################################################
# Simulation studies of time-to-event models with repeatability #
#################################################################

# load libraries
library(pacman) # checks if package is installed, if not installs it
p_load(survival, coxme, lme4, here, dplyr, purrr)

# load functions
source("R/function.R")

# from: https://stackoverflow.com/a/24569739
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}



##### load sim conditions -----------------------------------------------------

# load job array table
tab <- read.csv("sim/job_array.csv")

# get job number from pbs script
job <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
#job <- 5089

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




##### data generation  --------------------------------------------------------

# generate covariate (sex)
sex <- rbinom(n, size = 1, prob = 0.5)[rep(1:n_clusters, each = repetitions)]
# set true beta parameter value
beta <- 0.5

# generate cluster id variable
cluster <- rep(1:n_clusters, each = repetitions)

# generate frailty term (random effect) at the cluster level
frailty <- rnorm(n_clusters, sd = sqrt(sigma2.f))[cluster]

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

# get cut points for time intervals based on quantiles 
cutpoints_q <- quantile(dat$survival_time, probs = seq(0, 1, length.out = intervals + 1))
cutpoints_q <- cutpoints_q[-length(cutpoints_q)] # remove the last value

#### create a data frame (exploded_dat)
exploded_dat_q <- survSplit(Surv(survival_time, event) ~ sex + cluster,
                            data = dat,
                            cut = cutpoints_q, 
                            start = "tstart",
                            end = "tstop",
                            zero=0,
                            id = "observation")


# add time interval
exploded_dat_q$t_interval <- as.factor(exploded_dat_q$tstart)

# save current environment R objects
save.image(paste0("results/data/simdat_", seed, ".RDATA"))



##### modeling  -----------------------------------------------------------------------------------


###### coxme -------------
# fit a Cox prop hazards model with a frailty term using the coxme package

coxme_mod <- myTryCatch(coxme(Surv(survival_time, event) ~ sex + (1 | cluster), data = dat))

# save output
coxme_mod_res <- coxme_mod$value
coxme_mod_errors <- coxme_mod$error$message
coxme_mod_warnings <- coxme_mod$warning$message


# extract model estimates if there is no error 
if (is.null(coxme_mod_errors)) {
  
  # get model estimates
  beta_coxme <- fixef(coxme_mod_res)
  var_coxme <- VarCorr(coxme_mod_res)$cluster[[1]]
  ICC_coxme <- fr.lognormal(k, s, var_coxme, what = "tau")
  # get RMSE for each estimate
  beta_coxme_rmse <- beta - beta_coxme
  var_coxme_rmse <- sigma2.f - var_coxme
  # get error and warning messages to NA
  coxme_mod_errors <- NA
  coxme_mod_warnings <- ifelse(length(coxme_mod_warnings)==0, NA, paste(coxme_mod_warnings, collapse="; "))
  
# set to NA if there is an error message  
} else { 
  beta_coxme <- NA
  var_coxme <- NA
  ICC_coxme <- NA
  beta_coxme_rmse <- NA
  var_coxme_rmse <- NA
  coxme_mod_warnings <- ifelse(length(coxme_mod_warnings)==0, NA, paste(coxme_mod_warnings, collapse="; "))
}




#### coxph1 (normal) -------------
# fit a Cox prop hazards model with a frailty term using the survival package (assuming normal dist)

coxph1_mod <- myTryCatch(coxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gaussian"), dat))

# save output
coxph1_mod_res <- coxph1_mod$value
coxph1_mod_errors <- coxph1_mod$error$message
coxph1_mod_warnings <- coxph1_mod$warning$message


# extract model estimates if there is no error 
if (is.null(coxph1_mod_errors) & is.null(coxph1_mod_warnings)) {
  
  # get model estimates
  beta_coxph1 <- coxph1_mod_res$coefficients
  var_coxph_normal <- coxph1_mod_res$history$`frailty(cluster, distribution = "gaussian")`$theta
  ICC_coxph1 <- fr.lognormal(1, sigma2.f, var_coxph_normal, what = "tau")
  # get RMSE for each estimate
  beta_coxph1_rmse <- beta - beta_coxph1
  var_coxph1_rmse <- sigma2.f - var_coxph_normal
  # get error and warning messages to NA
  coxph1_mod_errors <- NA
  coxph1_mod_warnings <- ifelse(length(coxph1_mod_warnings)==0,NA, paste(coxph1_mod_warnings, collapse="; "))

# set to NA if there is an error
} else { 
  beta_coxph1 <- NA
  var_coxph_normal <- NA
  ICC_coxph1 <- NA
  beta_coxph1_rmse <- NA
  var_coxph1_rmse <- NA
  coxph1_mod_errors <- ifelse(length(coxph1_mod_errors)==0,NA, paste(coxph1_mod_errors, collapse="; "))
  coxph1_mod_warnings <- ifelse(length(coxph1_mod_warnings)==0,NA, paste(coxph1_mod_warnings, collapse="; "))
}




#### coxph2 (gamma) -------------
# fit a Cox prop hazards model with a frailty term using the survival package (assuming gamma dist)

coxph2_mod <-  myTryCatch(coxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gamma"), dat))

# save output
coxph2_mod_res <- coxph2_mod$value
coxph2_mod_errors <- coxph2_mod$error$message
coxph2_mod_warnings <- coxph2_mod$warning$message

# extract model estimates if there is no error 
if (is.null(coxph2_mod_errors)) {
  
  # get model estimates
  beta_coxph2 <- coxph2_mod_res$coefficients
  var_coxph_gamma <- coxph2_mod_res$history$`frailty(cluster, distribution = "gamma")`$theta
  ICC_coxph2 <- var_coxph_gamma/(var_coxph_gamma + 2)
  # get RMSE for each estimate
  beta_coxph2_rmse <- beta - beta_coxph2
  var_coxph2_rmse <- sigma2.f - var_coxph_gamma
  # get error and warning messages
  coxph2_mod_errors <- NA
  coxph2_mod_warnings <- ifelse(length(coxph2_mod_warnings)==0, NA, paste(coxph2_mod_warnings, collapse="; "))

# set to NA if there is an error message  
} else {
  beta_coxph2 <- NA
  var_coxph_gamma <- NA
  ICC_coxph2 <- NA
  beta_coxph2_rmse <- NA
  var_coxph2_rmse <- NA
  coxph2_mod_warnings <- ifelse(length(coxph2_mod_warnings)==0, NA, paste(coxph2_mod_warnings, collapse="; "))
}




#### glmm (binomial) -------------
# fit a discrete-time survival model with a frailty term using the glmer function from the lme4 package

glmm_mod <- myTryCatch(glmer(event ~ -1 + t_interval + sex + (1|cluster), 
                             data=exploded_dat_q, family=binomial(link="cloglog")))

# save output
glmm_mod_res <- glmm_mod$value
glmm_mod_errors <- glmm_mod$error$message
glmm_mod_warnings <- glmm_mod$warning$message


# extract model estimates if there is no error
if (is.null(glmm_mod_errors)) {
  
  # get model estimates
  beta_glmm <- fixef(glmm_mod_res)["sex"]
  var_glmm <- VarCorr(glmm_mod_res)$cluster[[1]]
  ICC_glmm <- var_glmm/(var_glmm + pi^2/6)
  # get RMSE for each estimate
  beta_glmm_rmse <- beta - beta_glmm
  var_glmm_rmse <- sigma2.f - var_glmm
  # get error and warning messages
  glmm_mod_errors <- NA
  glmm_mod_warnings <- ifelse(length(glmm_mod_warnings)==0, NA, paste(glmm_mod_warnings, collapse="; "))

# set to NA if there is an error message  
} else {
  beta_glmm <- NA
  var_glmm <- NA
  ICC_glmm <- NA
  beta_glmm_rmse <- NA
  var_glmm_rmse <- NA
  glmm_mod_warnings <- ifelse(length(glmm_mod_warnings)==0, NA, paste(glmm_mod_warnings, collapse="; "))
}




##### save results  ----------------------------------------------

res <- data.frame(name = rep(name, 4),
                  scenario = rep(scen, 4),
                  simseed_id = rep(seed, 4),
                  n = rep(n, 4),
                  n_clusters = rep(n_clusters, 4),
                  beta = rep(beta, 4),
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
                  sigma2.f_rmse = c(var_coxme_rmse,
                                    var_coxph1_rmse,
                                    var_coxph2_rmse,
                                    var_glmm_rmse),
                  beta_est = c(beta_coxme,
                               beta_coxph1,
                               beta_coxph2,
                               beta_glmm),
                  beta_rmse = c(beta_coxme_rmse,
                                beta_coxph1_rmse,
                                beta_coxph2_rmse,
                                beta_glmm_rmse),
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
