################################################################################
########################### Dubly Robust Estimator #############################
################################################################################

DREst <- function(a, x_a = NULL, x_b_nuis, x_c_nuis, y, delta_c, alpha.nuis, beta.nuis, gamma.nuis,
                  max.step, thres, alpha.start, fun_tx = NULL){
  dr.est <- dr.estimate(a, x_a, x_b_nuis, x_c_nuis, y, delta_c, alpha.nuis, beta.nuis, gamma.nuis,
                         max.step, thres, alpha.start, fun_tx = fun_tx)
  dr.point <- dr.est$point.esti
  convergence <- dr.est$convergence
  S_mat <- dr.est$S_mat

  return(list(dr.point = dr.point, convergence = convergence,
              S_mat = S_mat))
}



### MLE for porpensity score
fit_ps <- function(a, x_c_nuis, conf.level = 0.95) {
  # Fit the logistic model
  fit <- glm(a ~ x_c_nuis - 1, family = binomial(link = "logit"))

  # Extract estimates and standard errors
  coefs <- summary(fit)$coefficients
  est <- coefs[, "Estimate"]
  se  <- coefs[, "Std. Error"]

  # Compute z-value for the desired confidence level
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha/2)

  # Wald-type confidence intervals on the link scale
  lower <- est - z * se
  upper <- est + z * se

  # Package results
  return(list(point.est = est, se.est = se, conf.lower = lower,
              conf.upper = upper))
}
