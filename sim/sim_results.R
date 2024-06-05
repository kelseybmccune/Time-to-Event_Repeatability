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
    facet_wrap(~sigma2.f, labeller=labeller(sigma2.f =
                                              c("1" = "sigma2.f = 1",
                                                "2" = "sigma2.f = 2",
                                                "3" = "sigma2.f = 3"))) +
    theme_bw() 
  if (save) {
    filename <- sprintf("sim/output/%s_%s.png", name, variable_to_plot)
    ggsave(filename, plot = gg, width = 8, height = 6)
  }
  return(gg)
}






################################################################################

##### int2_cens0
# remove all model results if at least one model did not converge without warnings or errors
results_int2_cens0 <- results %>%
  group_by(simseed_id, scenario) %>%
  mutate(nowarnings = sum(is.na(model_warnings))) %>%
  filter(nowarnings == 4) %>%
  filter(censoring_prop==0, intervals==2)


# Save each plot 
plot_res(results_int2_cens0, "beta_est", name="beta_est_int2_cens0")
plot_res(results_int2_cens0, "sigma2.f_est", name="sigma2_est_int2_cens0")
plot_res(results_int2_cens0, "ICC", name="ICC_est_int2_cens0")


##### int4_cens0
# filter by censoring prop and interval split
results_int4_cens0 <- results %>%
  group_by(simseed_id, scenario) %>%
  mutate(nowarnings = sum(is.na(model_warnings))) %>%
  filter(nowarnings == 4) %>%
  filter(censoring_prop==0, intervals==4)


# Save each plot 
plot_res(results_int4_cens0, "beta_est", name="beta_est_int4_cens0")
plot_res(results_int4_cens0, "sigma2.f_est", name="sigma2_est_int4_cens0")
plot_res(results_int4_cens0, "ICC", name="ICC_est_int4_cens0")



##### int2_cens15
# filter by censoring prop and interval split
results_int2_cens15 <- results %>%
  group_by(simseed_id, scenario) %>%
  mutate(nowarnings = sum(is.na(model_warnings))) %>%
  filter(nowarnings == 4) %>%
  filter(censoring_prop==0.15, intervals==2)


# Save each plot
plot_res(results_int2_cens15, "beta_est", name="beta_est_int2_cens15")
plot_res(results_int2_cens15, "sigma2.f_est", name="sigma2_est_int2_cens15")
plot_res(results_int2_cens15, "ICC", name="ICC_est_int2_cens15")


##### int4_cens15
# filter by censoring prop and interval split
results_int4_cens15 <- results %>%
  group_by(simseed_id, scenario) %>%
  mutate(nowarnings = sum(is.na(model_warnings))) %>%
  filter(nowarnings == 4) %>%
  filter(censoring_prop==0.15, intervals==4)


# Save each plot 
plot_res(results_int4_cens15, "beta_est", name="beta_est_int4_cens15")
plot_res(results_int4_cens15, "sigma2.f_est", name="sigma2_est_int4_cens15")
plot_res(results_int4_cens15, "ICC", name="ICC_est_int4_cens15")






################################################################################

# get data for each scenario 1
results_scen1 <- results %>%
  group_by(simseed_id, scenario) %>%
  mutate(nowarnings = sum(is.na(model_warnings))) %>%
  filter(nowarnings == 4) %>%
  filter(scenario==1)

# Save each plot of performance measure
plot_res(results_scen1 , "beta_est", name="scen1")
plot_res(results_scen1 , "sigma2.f_est", name="scen1")
plot_res(results_scen1 , "ICC", name="scen1")

# get data for each scenario 2
results_scen2 <- results %>%
  group_by(simseed_id, scenario) %>%
  mutate(nowarnings = sum(is.na(model_warnings))) %>%
  filter(nowarnings == 4) %>%
  filter(scenario==2)

# Save each plot of performance measure
plot_res(results_scen2 , "beta_est", name="scen2")
plot_res(results_scen2 , "sigma2.f_est", name="scen2")
plot_res(results_scen2 , "ICC", name="scen2")


# get data for each scenario 3
results_scen3 <- results %>%
  group_by(simseed_id, scenario) %>%
  mutate(nowarnings = sum(is.na(model_warnings))) %>%
  filter(nowarnings == 4) %>%
  filter(scenario==3)

# Save each plot of performance measure
plot_res(results_scen3 , "beta_est", name="scen3")
plot_res(results_scen3 , "sigma2.f_est", name="scen3")
plot_res(results_scen3 , "ICC", name="scen3")

# get data for each scenario 4
results_scen4 <- results %>%
  group_by(simseed_id, scenario) %>%
  mutate(nowarnings = sum(is.na(model_warnings))) %>%
  filter(nowarnings == 4) %>%
  filter(scenario==4)

# Save each plot of performance measure
plot_res(results_scen4 , "beta_est", name="scen4")
plot_res(results_scen4 , "sigma2.f_est", name="scen4")
plot_res(results_scen4 , "ICC", name="scen4")



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

