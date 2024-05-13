# script for getting ICC from the coxme object

library(coxme)
library(here)
# example model

# we need an example data to create a function
dat <- read.csv(here("data","CTWemergence.csv"))
# "HT" variable indicates hiding time, or the latency to emerge; "Whorls" is a visual indicator of age
# 30 worms received 4 trials per day, across 4 days for a total of 16 trials. 

# 2 individuals have NA values, but it is not explained why. I'll assume these are censored (the worm didn't emerge in the trial time)
dat$event = ifelse(is.na(dat$HT),0,1)
dat$HT[which(is.na(dat$HT))]<- 375

dat$obs <- 1:nrow(dat)
#model <- coxme(Surv(HT, event)~ Whorls + (1|Worm_ID) + (1|obs), data=dat)
model <-  coxme(Surv(HT, event)~ Whorls + (1|Worm_ID), data=dat)

# Define the function
# missing values ignored

#' @title comxe_pval
#' @description This function calculates the p-value of the effect of the random effect in a coxme model. It also provides the p-value of the effect of the random effect using a bootstrapped method.
#' @param model A coxme model object
#' @param data The original data used to fit the model
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A vector of p-values
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author etc

coxme_pval <- function(model, data, boot = NULL) {
  # Get the original data
  
  if(all(class(model) %in% c("coxme")) == FALSE) {stop("Sorry, you need to fit a coxme model of class coxme")}
  
  # I think we need to use get the dimension of the data

  response <- as.data.frame(model$y[,1:2])
  
  fixed_formula <- as.formula(model$formulaList$fixed)
  
  # fit the model without any random effects
  fit <- survival::coxph(as.formula(fixed_formula), data = data)
  # loglikelihood ratio test
  # this is p value of effect of taking all random effects
  pval<- anova(fit, model)$P[-1]
  names(pval) <- "liklihood_ratio_test"
  
  if(!is.null(boot)){
  
  # we need to use replicate to create many vectors of these - randomize the data
  orders <- replicate(boot, sample(1:nrow(response)))  
  
  fixed_formula <- as.character(fixed_formula)  
  random_formula <-  as.vector(as.character(model$formulaList$random))        
  formula <- as.formula(paste("Surv(new_time, new_status)", 
                           "~", 
                           fixed_formula[3], 
                           "+",
                           paste(random_formula, collapse = "+")))
  data2 <- data

  # randomizaton/permutation tests
  pb <- progress::progress_bar$new(total = boot,
                                 format = "Bootstrapping [:bar] :percent ETA: :eta",
                                 show_after = 0)

  # loop
  num <- length(summary(model)$random$variance)

  store <- matrix(NA, nrow = num, ncol = boot)

  # Loop over the number of bootstraps
  for (i in 1:boot) {
  # Permute the data

  data2$new_time <- response$time[orders[ ,i]]
  data2$new_status <- response$status[orders[ ,i]]

  # Fit the original coxme model
  temp  <- tryCatch(coxme(formula, data = data2))

  # get variance component
   store[ ,i] <- summary(temp)$random$variance

   pb$tick()
   Sys.sleep(1 / boot)

    }
   
  # getting the p value
   pval2 <- sapply(1:num, function(x) {
     sum(store[x,] > summary(model)$random$variance[x])/boot}
     )
   
   names(pval2) <- paste(rep("bootstrapped_pval", num), 1:num, sep = "_")
  }

  if(exists("pval2")) {
    
  res <- c(pval, pval2)
  return(res)
  
  } else {
  res <- pval
  return(res)
  }

}


# test
coxme_pval(model, dat, boot = 1000)

#' @title coxme_icc_ci
#' @description This function calculates the 95 percent confidence interval for the intraclass correlation from the `coxme` objects.
#' @param model A coxme model object
#' @param upper.multiplier The multiplier for the upper bound of the confidence interval. Default is 10 (adjust to a higer value if the upper bound is not reached).
#' @return A vector of the lower, point estimate, and upper bounds of the 95 percent confidence interval for the intraclass correlation
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author etc

coxme_icc_ci <- function(model, upper.multiplier = 10) {
  if(all(class(model) %in% c("coxme")) == FALSE)
    {stop("Sorry, you need to fit a coxme model of class coxme")} 
  if(any(length(summary(model)$random$variance) > 1)) {stop("Sorry. At the moment, we can only have a model with one random effect.")}
  
  # Define a sequence of variance values``
  # the length of the response
  n <- nrow(model$y)
  cut = 100
  
  var_point <- summary(model)$random$variance
  
  # based on this pdf: https://cran.r-project.org/web/packages/coxme/vignettes/coxme.pdf
  # upper CI is limited to var_point*(10*log(n)) - so this could fail
  estvar1 <- seq(0.000001, var_point, length = cut)
  estvar2 <- seq(var_point, var_point*upper.multiplier, length = cut+1)[-1]
  estvar <- c(estvar1, estvar2)
  
  # Initialize a vector to store the log-likelihood values
  loglik <- double(cut)
  
  # Loop over the variance values
  for (i in seq_len(cut*2)) {
    # Fit a coxme model with fixed variance
    tfit <- update(model, vfixed = estvar[i])
    
    # Compute the log-likelihood
    loglik[i] <- 2 * diff(tfit$loglik)[1]
  }
  
  # Compute the threshold for the 95% confidence interval
  temp <-  as.numeric(2 * diff(model$loglik)[1]) - loglik
  
  # Find the variance values that correspond to the threshold
  # getting lower and upper CI using profile likelihood
  lower <- approx(temp[1:(cut)], sqrt(estvar[1:(cut)]), qchisq(.95, 1))$y
  upper <- approx(temp[(cut + 1):(2*cut)], sqrt(estvar[(cut + 1):(2*cut)]), qchisq(.95, 1))$y
  
  # Return the 95% confidence interval
  ICC_lower <- lower^2 / (lower^2 + pi^2 / 6)
  ICC_point <- var_point / (var_point + pi^2 / 6)
  ICC_upper <- upper^2 / (upper^2 + pi^2 / 6)
  names(ICC_lower) <- "lower"
  names(ICC_point) <- "ICC"
  names(ICC_upper) <- "upper"
  
  return(c(ICC_lower, ICC_point, ICC_upper))
}

# test
coxme_icc_ci(model)
