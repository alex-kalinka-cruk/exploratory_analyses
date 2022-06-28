# Helper functions.
.calc_beta1 <- function(lfc){
  # lfc - time-ordered log-fold change values for a single replicate.
  n <- length(lfc)
  x <- lfc
  xm <- c(0,lfc[1:(n-1)])
  prod_sum <- sum(x*xm)
  sum_prod <- sum(x)*sum(xm)
  sq_sum <- sum(xm^2)
  sum_sq <- sum(xm)^2
  num <- (1/n)*prod_sum - (1/n^2)*sum_prod
  denom <- (1/n)*sq_sum - (1/n^2)*sum_sq
  return(num/denom)
}

.calc_theta <- function(lfc, beta1, ts){
  # ts - time step absolute values in any units, e.g. 1, 3, 5 days.
  # Use mean time-step size.
  delta <- mean(diff(ts))
  theta <- -(1/delta)*log(beta1)
  return(theta)
}

.calc_mu <- function(lfc, beta1){
  n <- length(lfc)
  x <- lfc
  xm <- c(0,lfc[1:(n-1)])
  num <- (1/n)*sum(x - beta1*xm)
  denom <- 1 - beta1
  return(num/denom)
}

.calc_sigma <- function(lfc, beta1, theta, mu){
  n <- length(lfc)
  x <- lfc
  xm <- c(0,lfc[1:(n-1)])
  beta3 <- (1/n)*sum((x - beta1*xm - mu*(1 - beta1))^2)
  sigma <- sqrt(2*theta*beta3*(1/(1-beta1^2)))
  return(sigma)
}


#' 
#' 
#' 
calc_MLE_OU_lfc <- function(lfc, ts){
  beta1 <- .calc_beta1(lfc)
  theta <- .calc_theta(lfc, ts)
  mu <- .calc_mu(lfc, beta1)
  sigma <- .calc_sigma(lfc, beta1, theta, mu)
  return(data.frame(theta = theta, mu = mu, sigma = sigma))
}

