# Helper functions.
.calc_L_first_term <- function(x0, theta, mu, sigma){
  ft <- stats::dnorm((1/sigma)*sqrt(2*theta)*(x0-mu))
  return(ft)
}

.calc_L_second_term <- function(lfc, ts, theta, mu, sigma){
  n <- length(lfc)
  x <- lfc
  xm <- c(0,lfc[1:(n-1)])
  # Use mean time step delta.
  delta <- mean(diff(ts))
  mu_cond <- xm*exp(-theta*delta) + mu*(1-exp(-theta*delta))
  wb <- (1/sigma)*sqrt(2*theta*(1/(1-exp(-2*theta*delta))))*(x - mu_cond)
  st <- prod(stats::dnorm(wb))
  return(st)
}


#'
#'
calc_OU_likelihood <- function(lfc, ts, x0, theta, mu, sigma){
  st <- .calc_L_second_term(lfc, ts, theta, mu, sigma)
  return(st)
}
