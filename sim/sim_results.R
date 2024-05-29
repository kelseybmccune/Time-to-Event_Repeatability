##############################################
# Combine simulation results (on HPC/katana) #
##############################################

home.wd <- "/srv/scratch/z5394590/survival_repeatability/"

# initialise the result storage
results <- NULL

# change into appropriate result folder
setwd(paste(home.wd, "results/raw", sep = "/"))

# get a list of all the individual result files
res.list <- list.files()[grepl("res_", list.files(), fixed = T)]

# inner loop through individual sim-by-model files
for (job in res.list) {
  # get the job number
  job.no <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(job, ".", fixed = T), function(x){x[1]})), "_", fixed = T), function(x){x[length(x)]})))
  # load in the individual simulation results
  load(job)
  # rbind to results data frame
  results <- rbind(results, res)
  # remove the current results to ensure no duplicates
  rm(res)
}


# save the single result data frame within the results folder
setwd(home.wd)
save(list = "results", file = "collated_sim_results.RDATA")



##### Combine all RDATA file of dataframes into one list

# define the directory containing the Rdata files
folder_path <- "sim/output/data"

# get a list of all Rdata files in the folder
rdata_files <- list.files(folder_path, pattern = "simdat_\\d+\\.RDATA", full.names = T)

# initialise empty lists to store the data frames
dat_list <- list()
exploded_dat_q_list <- list()

# loop through each Rdata file and load the data frames and combine into lists
for (file in rdata_files) {
  load(file)
  if (exists("dat")) {
    dat_list[[file]] <- dat
  }
  if (exists("exploded_dat_q")) {
    exploded_dat_q_list[[file]] <- exploded_dat_q
  }
}

# Save the lists of data frames into a single Rdata file
save(dat_list, exploded_dat_q_list, file = "collated_sim_data.RDATA")




################
# Plot results #
################

load("sim/output/140524/collated_sim_results_subsets.RDATA")

# get job_number in scen.tab for results based on  scenario and simseed.id in the results dataset
scen <- scen.tab %>%
  select(scenario,
         simseed_id=sim,
         job_number)

t <- merge(results, scen, by=c("simseed_id", "scenario"))


# check which sequences of sim datasets have missing numbers
sequence <- 1:12000
setdiff(sequence, t$job_number)

# 5089 11089


# load libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggdark)
library(ggplot2)
library(ggdist)  # for half-eye plots
library(gridExtra)
library(ggsave)




################################################################################

# Function to plot performance measures ------------

plot_res <- function(res, variable_to_plot, name="res", save=TRUE) {
  
  # reorder models for plots
  res$model <- factor(res$model, levels=c("coxme",
                                          "coxph_normal",
                                          "glmm",
                                          "coxph_gamma"))
  
  # create half eye plot
  gg <- ggplot(res, aes(x=model, y=get(variable_to_plot), color=model, fill=model)) +
    stat_halfeye() +  
    geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=0.6) + 
    geom_boxplot(width=0.1, alpha=0.4) +
    scale_color_manual(values=c("#CE72DD", "#FECA91", "#73B496", "#8C7CBB")) +
    scale_fill_manual(values=alpha(c("#CE72DD", "#FECA91", "#73B496", "#8C7CBB"), 0.4)) + 
    labs(title=sprintf("Plot of %s", variable_to_plot),x="model",y=variable_to_plot) +
    facet_wrap(~sigma2.f, labeller=labeller(sigma2.f =
                                              c("1" = "sigma2.f = 1",
                                                "2" = "sigma2.f = 2",
                                                "3" = "sigma2.f = 3"))) +
    theme_bw() 
  if (save) {
    filename <- sprintf("sim/output/%s_%s.png", name, variable_to_plot)
    ggsave(filename, plot = gg, width = 8, height = 5)
  }
  return(gg)
}


################################################################################

# filter by censoring prop and interval split
results_int2_cens0 <- results %>%
  filter(censoring_prop==0, intervals==2)


# filter all simseed.id without at least one model_warnings
results_int2_cens0 <- results_int2_cens0 %>%
  group_by(simseed_id, sigma2.f) %>%
  mutate(nowarnings = sum(is.na(model_warnings))) %>%
  filter(nowarnings == 4)


# Save each plot of performance measure using par(mfrow=c(1,3))
par(mfrow=c(1,3))
plot_res(results_int2_cens0, "beta_est", name="beta_est_int2_cens0")
plot_res(results_int2_cens0, "sigma2.f_est", name="sigma2_est_int2_cens0")
plot_res(results_int2_cens0, "ICC", name="ICC_est_int2_cens0")


# filter by censoring prop and interval split
results_int4_cens0 <- results %>%
  filter(censoring_prop==0, intervals==4)


# Save each plot of performance measure using par(mfrow=c(1,3))
plot_res(results_int4_cens0, "beta_est", name="beta_est_int4_cens0")
plot_res(results_int4_cens0, "sigma2.f_est", name="sigma2_est_int4_cens0")
plot_res(results_int4_cens0, "ICC", name="ICC_est_int4_cens0")







################################################################################

# get data for each scenario 1
results_scen1 <- results %>%
  filter(scenario==1)
  
# Save each plot of performance measure
plot_res(results_scen1 , "beta_est", name="scen1")
plot_res(results_scen1 , "sigma2.f_est", name="scen1")
plot_res(results_scen1 , "ICC", name="scen1")

# get data for each scenario 2
results_scen2 <- results %>%
  filter(scenario==2)

# Save each plot of performance measure
plot_res(results_scen2 , "beta_est", name="scen2")
plot_res(results_scen2 , "sigma2.f_est", name="scen2")
plot_res(results_scen2 , "ICC", name="scen2")










################################################################################


# # filter out simseed.id if there is at least one warning from any model
# results_scen1_red <- results_scen1 %>%
#   group_by(simseed_id) %>%
#   mutate(nowarnings = sum(is.na(model_warnings))) %>%
#   filter(nowarnings == 4)
# 
# # Save each plot of performance measure
# plot_res(results_scen1_red, "beta_est", name="nowarnings")
# plot_res(results_scen1_red, "sigma2.f_est", name="nowarnings")
# plot_res(results_scen1_red, "ICC", name="nowarnings")

