# simulation

# Load the necessary packages
library(survival)
library(coxme)
library(lme4)
library(here)

#loading function.R 

source(here("R","function.R"))

set.seed(777) # for reproducibility

# Number of observations
n <- 500 # the number of observations
n_clusters <- 100 # the number of individuals (cluster)

# n = 30 n_cluster = 15


repetitions <- n/n_clusters

# Generate some covariates
#age <- rnorm(n, mean = 50, sd = 10)
# we need to assign sex to cluster (individual not observation)
sex <- rbinom(n, size = 1, prob = 0.5)
# True parameter values
beta <- 2

# Generate cluster variable (e.g., 10 clusters)
cluster <- rep(1:n_clusters, each = repetitions)
#cluster <- sample(1:n_clusters, n, replace = TRUE)

# Generate frailty term at the cluster level
# 4 senarios var = 0.5, 2, 4
# this is the random effect
frailty <- rnorm(n_clusters, sd = sqrt(2))
frailty <- frailty[cluster]

# Generate survival times from an exponential distribution
hazard <- exp(beta * sex + frailty)
survival_time <- rexp(n, hazard)

# # we could use Weibull - probably not 
#shape <- 1  # shape parameter for the Weibull distribution
#scale <- 1/hazard  # scale parameter for the Weibull distribution
#survival_time <- rweibull(n, shape, scale)

# Introduce censoring for 15% of the individuals
# we could model this differently
# censoring 0 or 15??
censoring_prop <- 0.15 # 15%
censoring_time <- quantile(survival_time, 1 - censoring_prop)


event <- ifelse(survival_time <= censoring_time, 1, 0)
survival_time <- pmin(survival_time, censoring_time)

# Create a data frame
dat <- data.frame(observation = 1:n, cluster = cluster, sex = sex, survival_time = survival_time, event = event)

# Print the summary of the model
#summary(coxme_fit)
#VarCorr(coxme_fit)$cluster[[1]]

# this can be 3 or 9     
intervals <- 9 # the number of intervals

# Fit a discrete-time survival model with a frailty term using the glmer function from the lme4 package
exploded_dat <- survSplit(Surv(survival_time, event) ~ sex + cluster, data = dat,
                          cut = (1:intervals)*(max(dat$survival_time)/intervals), 
                          start = "tstart",end = "tstop",  zero=0, id = "observation")

# Poisson offset - not needed
#exploded_dat$os_time <- exploded_dat$tstop - exploded_dat$tstart

# time interval
exploded_dat$t_interval <- as.factor(exploded_dat$tstart)


# modeling
# Fit a Cox proportional hazards model with a frailty term using the coxme package
coxme_fit <- coxme(Surv(survival_time, event) ~ sex + (1|cluster), data = dat)
coxph_fit1 <- coxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gaussian"), dat)
coxph_fit2 <- coxph(Surv(survival_time, event) ~ sex + frailty(cluster, distribution="gamma"), dat)
glmm_fit <- glmer(event~ -1 + t_interval + sex + (1|cluster), 
              data=exploded_dat, family=binomial(link="cloglog"))

var_coxme <- VarCorr(coxme_fit)$cluster[[1]]
var_coxph_normal <- coxph_fit1$history$`frailty(cluster, distribution = "gaussian")`$theta
var_glmm <- VarCorr(glmm_fit)$cluster[[1]]
var_coxph_gamma <- coxph_fit2$history$`frailty(cluster, distribution = "gamma")`$theta # expected to be different

var_coxme
var_coxph_normal
var_glmm
var_coxph_gamma
# fixed effects
fixef(coxme_fit)
fixef(glmm_fit)["sex"]
coxph_fit1$coefficients
coxph_fit2$coefficients

# ICC tau and glmmer
fr.lognormal(k,s,var_coxme, what = "tau")
fr.lognormal(k,s,var_coxph_normal, what = "tau")
var_glmm/(var_glmm + pi^2/6) # expected to be different
var_coxph_gamma/(var_coxph_gamma + 2)
