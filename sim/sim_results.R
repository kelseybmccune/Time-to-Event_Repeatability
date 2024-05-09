##########################################
# Combine simulation results (on katana) #
##########################################

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
  results <- rbind(dat, res)
  # remove the current results to ensure no duplicates
  rm(res)
}


# save the single result data frame within the results folder
setwd(home.wd)
save(list = "results", file = "collated_sim_results.RDATA")



################
# Plot results #
################

load("sim/output/collated_sim_results.RDATA")

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
    theme_bw() 
  if (save) {
    filename <- sprintf("sim/output/%s_%s.png", name, variable_to_plot)
    ggsave(filename, plot = gg, width = 6, height = 4)
  }
  return(gg)
}


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



###############################################################################

##### split by scenario number


# Plot of beta_est vs models
beta_est_plot <- ggplot(results, aes(x=model, y=beta_est, color=model, fill=model)) + 
  stat_halfeye() +
  scale_color_manual(values=c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"))+
  scale_fill_manual(values=alpha(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 0.4)) +
  labs(title="Beta estimate",x="Model", y = "Beta estimate")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=2, colour="darkgray")+ # beta=2
  theme_bw()
beta_est_plot


# Plot of sigma2.frailty vs models
sigma2_est_plot <- ggplot(results, aes(x=model, y=sigma2.f_est, color=model, fill=model)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  labs(title="Variance estimate",x="Model", y = "Variance estimate")+
  geom_boxplot(width=0.4)+
  theme_bw()
sigma2_est_plot


# Plot of ICC vs models
sigma2_est_plot <- ggplot(results, aes(x=model, y=ICC, color=model, fill=model)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  labs(title="ICC estimate",x="Model", y = "ICC estimate")+
  geom_boxplot(width=0.4)+
  theme_bw()
sigma2_est_plot

