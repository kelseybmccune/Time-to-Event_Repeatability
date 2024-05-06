# function
# obtained from
# https://github.com/cran/parfm/blob/master/R/frailtydistros.R

g <- function(w, k, s, sigma2) {
  -k * w + exp(w) * s + w ^ 2 /  (2 * sigma2)
}

g1 <- function(w, k, s, sigma2) {
  -k + exp(w) * s + w / sigma2
}

g2 <- function(w, k, s, sigma2) {
  exp(w) * s + 1 / sigma2    
} 

Lapl <- Vectorize(function(s, k, sigma2) {
  # Find wTilde = max(g(w)) so that g'(wTilde; k, s, theta) = 0
  WARN <- getOption("warn")
  options(warn = -1)
  wTilde <- optimize(f = g, c(-1e10, 1e10), maximum = FALSE,
                     k = k, s = s, sigma2 = sigma2)$minimum
  options(warn = WARN)
  
  # Approximate the integral via Laplacian method
  res <- (-1) ^ k * 
    exp(-g(w = wTilde, k = k, s = s, sigma2 = sigma2)) /
    sqrt(sigma2 * g2(w = wTilde, k = k, s = s, sigma2 = sigma2))
  return(res)
}, 's')

intTau <- Vectorize(function(x, intTau.sigma2=sigma2) {
  res <- x * 
    Lapl(s = x, k = 0, sigma2 = intTau.sigma2) *
    Lapl(s = x, k = 2, sigma2 = intTau.sigma2)
  return(res)
}, "x")

fr.lognormal <- function(k,
                         s,
                         sigma2,
                         what = "logLT") {
  # if (!(is.numeric(sigma2) && (sigma2 > 0)))
  # stop("The parameter sigma2 is not a positive value.")
  
  if (what == "logLT") {
    # if (!(is.numeric(s) && (s > 0)))
    #     stop("The parameter s is not positive.")
    # Find wTilde = max(g(w)) so that g'(wTilde; k, s, theta) = 0
    WARN <- getOption("warn")
    options(warn = -1)
    wTilde <- nlm(f = g, p = 0, k = k, s = s, sigma2 = sigma2)$estimate
    options(warn = WARN)
    
    # Approximate the integral via Laplacian method
    res <- -g(w = wTilde, k = k, s = s, sigma2 = sigma2) -
      log(sigma2 * g2(w = wTilde, k = k, s = s, sigma2 = sigma2)
      ) / 2
    return(res)
  }
  else if (what == "tau") {
    intTau <- Vectorize(function(x, intTau.sigma2=sigma2) {
      res <- x * 
        Lapl(s = x, k = 0, sigma2 = intTau.sigma2) *
        Lapl(s = x, k = 2, sigma2 = intTau.sigma2)
      return(res)
    }, "x")
    
    tauRes <- 4 * integrate(
      f = intTau, lower = 0, upper = Inf, 
      intTau.sigma2 = sigma2)$value - 1
    return(tauRes)
  }
}
