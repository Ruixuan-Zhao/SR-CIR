################################################################################
####################### Maximum likelihood: CIR + SOP ##########################
################################################################################

cir_loglik_i <- function(pars, a, fun_tx, x_sop, y, delta_c, timepoints) {
  n <- length(a)
  m <- length(timepoints)
  p_rho <- fun_tx_dim.cir(fun_tx)
  p_tau <- ncol(x_sop)
  rho <- pars[seq_len(p_rho)]
  tau <- pars[p_rho + seq_len(p_tau)]
  p0p1 <- getProb.sequent.cir(rho, tau, fun_tx, x_sop, timepoints)

  ll <- numeric(n)
  for (i in seq_len(n)) {
    k <- match(y[i], timepoints)
    arm <- a[i] + 1L
    tail_probs <- vapply(k:m, function(j) p0p1[[j]][i, arm], numeric(1))
    F_curr <- prod(tail_probs)
    if (k == 1L) {
      mass <- F_curr
    } else {
      F_prev <- F_curr * p0p1[[k - 1L]][i, arm]
      mass <- F_curr - F_prev
    }
    surv_after <- 1 - F_curr
    if (delta_c[i] == 1) {
      ll[i] <- log(mass)
    } else {
      ll[i] <- log(surv_after)
    }
  }
  ll
}

max.likelihood.cir <- function(a, fun_tx, x_sop, y, delta_c, rho.start,
                               tau.start, max.step = 1000, thres = 1e-6,
                               timepoints = NULL) {
  x_sop <- as_matrix.cir(x_sop, name = "x_sop")
  n <- length(a)
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  p_rho <- fun_tx_dim.cir(fun_tx)
  p_tau <- ncol(x_sop)
  if (is.null(rho.start)) rho.start <- rep(0, p_rho)
  if (is.null(tau.start)) tau.start <- rep(0, p_tau)

  neg_loglik <- function(pars) {
    -sum(cir_loglik_i(pars, a, fun_tx, x_sop, y, delta_c, timepoints))
  }

  rho <- rho.start
  tau <- tau.start
  diff <- thres + 1
  step <- 0L
  opt_rho <- opt_tau <- list(convergence = 0L)

  rel_diff <- function(new, old) {
    sum((new - old)^2) / sum(new^2 + thres)
  }
  neg_loglik_rho <- function(rho_try) {
    neg_loglik(c(rho_try, tau))
  }
  neg_loglik_tau <- function(tau_try) {
    neg_loglik(c(rho, tau_try))
  }

  while (diff > thres && step < max.step) {
    step <- step + 1L

    opt_rho <- stats::optim(
      rho, neg_loglik_rho,
      method = "BFGS",
      control = list(maxit = max(100, floor(max.step / 10)), reltol = thres)
    )
    diff_rho <- rel_diff(opt_rho$par, rho)
    rho <- opt_rho$par

    opt_tau <- stats::optim(
      tau, neg_loglik_tau,
      method = "BFGS",
      control = list(maxit = max(100, floor(max.step / 10)), reltol = thres)
    )
    diff_tau <- rel_diff(opt_tau$par, tau)
    tau <- opt_tau$par

    diff <- max(diff_rho, diff_tau)
  }

  point <- c(rho, tau)
  hessian <- numDeriv::hessian(neg_loglik, point)
  score_i <- numDeriv::jacobian(function(pars) {
    cir_loglik_i(pars, a, fun_tx, x_sop, y, delta_c, timepoints)
  }, point)
  cov <- safe_solve.cir(hessian)
  convergence <- if (step < max.step &&
                     identical(opt_rho$convergence, 0L) &&
                     identical(opt_tau$convergence, 0L)) 0L else 1L
  list(
    par = point,
    convergence = convergence,
    value = neg_loglik(point),
    cov = cov,
    S_mat = score_i,
    step = step,
    timepoints = timepoints
  )
}
