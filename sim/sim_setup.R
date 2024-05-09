##############################################
# Simulation conditions and job array set up #
##############################################

name <- "test"

#### model parameters
n <- c(1000)                    # number of observations
n_clusters <- c(200, 500)       # number of individuals (cluster)
sigma2.f <- c(1, 4)             # variance of frailty random effect (removed 0.5 for now)
censoring_prop <-  c(0.15)  # proportion of censored data
intervals <- c(2,4)           # number of time intervals

# number of replications per scenario
repl <- 100
sim <- rep(1:repl)


#### scenario and job array table set up 

# make table of all scenarios
scen.tab <- expand.grid(sim = sim, name = name, 
                        n = n, n_clusters = n_clusters, 
                        sigma2.f = sigma2.f, 
                        censoring_prop = censoring_prop,
                        intervals = intervals)

# # keep specific combinations of n and n_clusters
# scen.tab <- scen.tab %>%
#   filter((n == 30 & n_clusters == 15) | (n == 500 & n_clusters == 100))

# add columns for to store results and job number
scen.tab$save_location <- rep("/srv/scratch/z5394590/survival_repeatability/", each=nrow(scen.tab))
scen.tab$job_number <- c(1:nrow(scen.tab))
conditions <- 2 * length(sigma2.f) * length(censoring_prop) * length(intervals)
scen.tab$scenario <- rep(1:conditions, each = repl)

# save as csv file
write.csv(scen.tab, "sim/job_array.csv", row.names = FALSE)
