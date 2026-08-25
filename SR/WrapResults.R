################################################################################
######## Function for wrapping estimation results into a nice format ###########
################################################################################

# refer to R package brm

WrapResults = function(point.est, cov, S_mat, name, x_a, x_b, converged,
                       fun_tx = NULL) {

  se.est <- sqrt(diag(cov %*% t(S_mat) %*% S_mat %*% t(cov)))


  conf.lower = point.est + stats::qnorm(0.025) * se.est
  conf.upper = point.est + stats::qnorm(0.975) * se.est
  p.temp = stats::pnorm(point.est/se.est, 0, 1)
  p.value = 2 * pmin(p.temp, 1 - p.temp)

  names(point.est) = names(se.est) = rownames(cov) = colnames(cov) = names(conf.lower) = names(conf.upper) = names(p.value) = name

  coefficients = cbind(point.est, se.est, conf.lower, conf.upper, p.value)

  n_alpha = sum(grepl("^alpha", name))
  alpha.est = point.est[seq_len(n_alpha)]

  if (!is.null(fun_tx)) {
    param.est = do.call(cbind, lapply(fun_tx, function(F_tx) {
      F_tx = if (is.null(dim(F_tx))) matrix(F_tx, ncol = 1) else as.matrix(F_tx)
      exp(F_tx %*% alpha.est)
    }))
  } else if (is.vector(x_a)){
    linear.predictors = x_a * point.est[1]
    param.est = exp(linear.predictors)
  }else{
    linear.predictors = x_a %*% point.est[1:ncol(x_a)]
    param.est = exp(linear.predictors)
  }

  sol = list(point.est = point.est, se.est = se.est, cov = cov, S_mat = S_mat,
             conf.lower = conf.lower, conf.upper = conf.upper, p.value = p.value,
             coefficients = coefficients, param.est = param.est, x_a = x_a, x_b = x_b,
             fun_tx = fun_tx, converged = converged)
  class(sol) = c("brm", "list")
  attr(sol, "hidden") = c("se.est", "cov", "conf.lower", "conf.upper",
                          "p.value","coefficients", "param.est", "x_a", "x_b",
                          "fun_tx", "converged")

  return(sol)

}
