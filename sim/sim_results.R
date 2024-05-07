##############################
# Combine simulation results #
##############################

home.wd <- "/srv/scratch/z5394590/survival_repeatability/"

# initialise the result storage
dat <- NULL

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
  dat <- rbind(dat, res)
  # remove the current results to ensure no duplicates
  rm(res)
}


# save the single result data frame within the results folder
setwd(home.wd)
save(list = "dat", file = "collated_sim_results.RDATA")



################
# Plot results #
################

load("data/collated_sim_results.RDATA")

results <- dat
results$model <- factor(results$model, levels=c("coxme",
                                                "coxph_normal",
                                                "coxph_gamma",
                                                "glmm"))

# remove outlier points
results <- results %>%
  filter(beta_est < 5,
         sigma2.f_est < 10)

# load libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggdark)
library(ggplot2)
library(ggdist)  # for half-eye plots
library(gridExtra)
library(ggsave)

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


################################################################################

# Function to plot performance measures ------------

plot_results <- function(res, variable_to_plot, name="res", save=TRUE) {
  
  # Reorder models for plots
  res$model <- factor(res$model, levels=c("coxme",
                                          "coxph_normal",
                                          "glmm",
                                          "coxph_gamma"))
  
  # Create ggplot object
  gg <- ggplot(res, aes(x=model, y=get(variable_to_plot), color=model, fill=model)) +
    stat_halfeye() +  # Replace with half-eye plot from ggdist
    geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=0.6) +  # Add points
    geom_hline(aes(yintercept=0), color="white") + # Add line at zero
    geom_boxplot(width=0.1, alpha=0.4) +
    scale_color_manual(values=c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B")) +
    scale_fill_manual(values=alpha(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 0.4)) + # Fill with semi-transparent pastels
    labs(title=sprintf("Plot of %s", variable_to_plot),
         x="Package",
         y=variable_to_plot) +
    scale_shape_manual(name="Species Size", values=c(16, 17, 18, 19)) +
    facet_wrap(~ factor(species_size),
               labeller = labeller(species_size = function(x) paste("Species size:", x)))
  
  if (save) {
    filename <- sprintf("output/sim_simple/%s_%s_%diter_sim_simple.png",
                        name, variable_to_plot, iters)
    ggsave(filename, plot = gg, width = 8, height = 4)
  }
  
  return(gg)
}


# Save each plot of performance measure
plot_results(results, "run_time")
plot_results(results, "mu_bias")
plot_results(results, "mu_mse")
plot_results(results, "mu_ci_width")
plot_results(results, "s2_sp_bias")
plot_results(results, "s2_phylo_bias")
plot_results(results, "s2_resid_bias")
