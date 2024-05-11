# figure 2 - the relationship between ICCnp and ICC 

# packages

library(coxme)
library(here)
library(tidyverse)
#library(patchwork)

# load the function
source(here("R","function.R"))


# create a vector (sequence) of values form 0 to 5 divided by 0.1
sigma2 <- seq(0, 15, by = 0.05)[-1]

# use this sigma2 to calculate the ICCnp
ICCnp <- sapply(sigma2, function(x)fr.lognormal(k,s,x,what = "tau"))

ICC <- sapply(sigma2, function(x) x/(x + pi^2/6))
  
# create a data frame

df <- data.frame(sigma2, ICCnp, ICC)

# Reshape the data to long format
df_long <- tidyr::pivot_longer(df, c(ICCnp, ICC), names_to = "type", values_to = "value")

# Plot the data with a legend box inside the plot
ggplot(df_long, aes(x = sigma2, y = value, color = type)) +
  geom_line() +
  labs(#x = paste("Variance","(","sigma^2",")"), 
       x = expression(paste("Variance (", sigma[alpha]^2, ")")) ,
       y = "Repeatablity (intra-class correlation)",parse = TRUE) +
  scale_color_manual(values = c("blue", "red"), labels = c("ICC", "ICCnp")) +  
  ylim(0,1) +# Switched the labels
  theme_bw() +
  theme(legend.position = c(0.85, 0.15))+
  guides(color = guide_legend(title = NULL)) +
  theme(
    legend.background = element_rect(
      colour = "black",  # Change the color of the box
      fill = "white",  # Change the fill color of the box
      size = 0.2,  # Change the size of the box border
      linetype = "solid"  # Change the type of the line
    )
  )


# plot the data showing a relationship 
#between ICCnp and ICC not with sigma2 - also add 1:1 line (dotted line)
# ICC as y and ICCnp as x
# ggplot(df, aes(x = ICCnp, y = ICC)) +
#   geom_line(color = "blue") +
#   geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
#   labs(title = "The relationship between ICCnp and ICC",
#        x = "ICCnp",
#        y = "ICC") +
#   theme_minimal()



