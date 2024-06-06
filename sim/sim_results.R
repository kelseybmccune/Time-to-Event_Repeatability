################
# Plot results #
################

# load libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggdark)
library(ggplot2)
library(ggdist)  # for half-eye plots
library(gridExtra)
library(ggsave)
library(patchwork)


# get results from all sims
load("sim/output/collated_sim_results.RDATA")

# # get job_number in scen.tab for results based on  scenario and simseed.id in the results dataset
# scen <- scen.tab %>%
#   select(scenario,
#          simseed_id=sim,
#          job_number)
# 
# t <- merge(results, scen, by=c("simseed_id", "scenario"))
# 
# # check which sequences of sim datasets have missing numbers
# sequence <- 1:12000
# setdiff(sequence, t$job_number)
# # 5089 11089 iterations did not run 



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
    #facet_wrap(~sigma2.f, labeller=labeller(sigma2.f =
    #                                          c("1" = "sigma2.f = 1",
    #                                            "2" = "sigma2.f = 2",
    #                                            "3" = "sigma2.f = 3"))) +
    theme_bw()+
    theme(legend.position="none")
  
  # Add number of observations per model
  gg <- gg + geom_text(data = res %>% group_by(model) %>% summarise(count = n()), 
                       aes(label = paste0("n = ", count), y = Inf), vjust = 1.5, color = "black")
    
  # if (save) {
  #   filename <- sprintf("sim/output/plots/%s_%s.png", name, variable_to_plot)
  #   ggsave(filename, plot = gg, width = 7, height = 5)
  # }
  return(gg)
}


################################################################################

# Function to get data for each scenario ---------------
get_scenario_data <- function(results, scenario_num) {
  results %>%
    group_by(simseed_id, scenario) %>%
    mutate(nowarnings = sum(is.na(model_warnings))) %>% ## remove observations if at least one model had a warning 
    filter(nowarnings == 4) %>%
    filter(scenario == scenario_num)
}


################################################################################

# get labels for each scenario and variable (for plot title)
scenario_labels <- c("Interval = 2, Censoring = 0, sigma2.f = 1",
                     "Interval = 2, Censoring = 0, sigma2.f = 2",
                     "Interval = 2, Censoring = 0, sigma2.f = 3",
                     "Interval = 2, Censoring = 0.15, sigma2.f = 1",
                     "Interval = 2, Censoring = 0.15, sigma2.f = 2",
                     "Interval = 2, Censoring = 0.15, sigma2.f = 3",
                     "Interval = 4, Censoring = 0, sigma2.f = 1",
                     "Interval = 4, Censoring = 0, sigma2.f = 2",
                     "Interval = 4, Censoring = 0, sigma2.f = 3",
                     "Interval = 4, Censoring = 0.15, sigma2.f = 1",
                     "Interval = 4, Censoring = 0.15, sigma2.f = 2",
                     "Interval = 4, Censoring = 0.15, sigma2.f = 3")


# generate plot for each scenario and performance measure
for (scenario_num in 1:12) {
  results_scen <- get_scenario_data(results, scenario_num)
  
  beta_plot <- plot_res(results_scen, "beta_est")
  sigma2_plot <- plot_res(results_scen, "sigma2.f_est")
  ICC_plot <- plot_res(results_scen, "ICC")
  
  final_plot <- (beta_plot | sigma2_plot | ICC_plot) + 
    plot_layout(ncol = 3) + 
    plot_annotation(
      title = sprintf("%s (Scenario %d)", scenario_labels[scenario_num], scenario_num),
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
  
  filename <- sprintf("sim/output/plots/plot_scenario_%d.png", scenario_num)
  ggsave(filename, plot = final_plot, width = 12, height = 6)
  
  print(final_plot)
}





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

