################################################################################
############################ CIR MLE wrapper ###################################
################################################################################

MLEst.cir <- function(a, fun_tx, x_sop, y, delta_c, rho.start = NULL,
                      tau.start = NULL, max.step = 1000, thres = 1e-6,
                      timepoints = NULL) {
  x_sop <- as_matrix.cir(x_sop, name = "x_sop")
  n <- length(a)
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  check_fun_tx.cir(fun_tx, n = n, m = length(timepoints))
  p_rho <- fun_tx_dim.cir(fun_tx)
  p_tau <- ncol(x_sop)
  if (is.null(rho.start)) rho.start <- rep(0, p_rho)
  if (is.null(tau.start)) tau.start <- rep(0, p_tau)

  mle <- max.likelihood.cir(
    a = a, fun_tx = fun_tx, x_sop = x_sop, y = y, delta_c = delta_c,
    rho.start = rho.start, tau.start = tau.start,
    max.step = max.step, thres = thres, timepoints = timepoints
  )
  point.est <- mle$par
   cov.sandwich <- mle$cov %*% crossprod(mle$S_mat) %*% t(mle$cov)
  se.est <- sqrt(pmax(diag(cov.sandwich), 0))
  conf.lower <- point.est + stats::qnorm(0.025) * se.est
  conf.upper <- point.est + stats::qnorm(0.975) * se.est
  p.temp <- stats::pnorm(point.est / se.est, 0, 1)
  p.value <- 2 * ifelse(p.temp < 0.5, p.temp, 1 - p.temp)
  name <- c(paste0("rho", seq_len(p_rho) - 1L), paste0("tau", seq_len(p_tau) - 1L))
  names(point.est) <- names(se.est) <- names(conf.lower) <- names(conf.upper) <- names(p.value) <- name
  rownames(mle$cov) <- colnames(mle$cov) <- name
  list(
    point.est = point.est,
    se.est = se.est,
    cov = mle$cov,
    S_mat = mle$S_mat,
    conf.lower = conf.lower,
    conf.upper = conf.upper,
    p.value = p.value,
    convergence = mle$convergence,
    timepoints = mle$timepoints
  )
}
