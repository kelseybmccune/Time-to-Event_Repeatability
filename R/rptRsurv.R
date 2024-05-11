# script for getting ICC from the coxme object

library(coxme)


# Define the function
permute_coxme <- function(model) {
  # Get the original data
  data <- model$data
  
  # Randomize the response variable
  response <- data[, c("time", "event")]
  response <- response[sample(nrow(response)), ]
  
  # Update the data with the randomized response
  data$time <- response$time
  data$event <- response$event
  
  # Fit a new coxme model with the randomized data
  permuted_model <- coxme(Surv(time, event) ~ ., data = data, random = model$random)
  
  # Return the permuted model
  return(permuted_model)
}


# Fit the original coxme model
model <- coxme(Surv(time, event) ~ ., data = my_data, random = ~ 1 | group)

# Fit a coxme model with permuted data
permuted_model <- permute_coxme(model)

# Print the summary of the permuted model
summary(permuted_model)


# Install the necessary packages if not already installed
if (!require(coxme)) {
  install.packages("coxme")
}

# Load the necessary package
library(coxme)

# Define the function
get_ci_coxme <- function(model, n = 100) {
  # Define a sequence of variance values
  estvar <- seq(0.01, 1, length = n)^2
  
  # Initialize a vector to store the log-likelihood values
  loglik <- double(n)
  
  # Loop over the variance values
  for (i in seq_len(n)) {
    # Fit a coxme model with fixed variance
    tfit <- update(model, vfixed = estvar[i])
    
    # Compute the log-likelihood
    loglik[i] <- 2 * diff(tfit$loglik)[1]
  }
  
  # Compute the threshold for the 95% confidence interval
  threshold <- 2 * diff(model$loglik)[1] - qchisq(0.95, 1)
  
  # Find the variance values that correspond to the threshold
  lower <- approx(loglik[1:(n/2)], sqrt(estvar[1:(n/2)]), threshold)$y
  upper <- approx(loglik[(n/2 + 1):n], sqrt(estvar[(n/2 + 1):n]), threshold)$y
  
  # Return the 95% confidence interval
  return(c(lower, upper))
}

# Fit the original coxme model
model <- coxme(Surv(time, event) ~ . + (1 | group), data = my_data)

# Compute the 95% confidence interval for the standard deviation of the random effect
ci <- get_ci_coxme(model)

# Print the confidence interval
print(ci)