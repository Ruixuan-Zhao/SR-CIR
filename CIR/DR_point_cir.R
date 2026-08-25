################################################################################
######################## CIR doubly robust estimator ###########################
################################################################################

cir_build_Q.cir <- function(Fa) {
  m <- length(Fa)
  F_prev <- c(0, Fa[-m])
  F_next <- c(Fa[-1], 1)
  q1 <- -Fa * (1 - F_prev) / (Fa - F_prev)
  q2 <- Fa * (1 - F_next) / (F_next - Fa)
  q2[m] <- 0
  Gamma <- matrix(0, m, m)
  if (m > 1L) {
    for (k in seq_len(m - 1L)) Gamma[k, k + 1L] <- 1
  }
  diag(q1, m, m) + diag(q2, m, m) %*% Gamma
}

cir_R.cir <- function(Fa, y_i, delta_i, timepoints) {
  m <- length(timepoints)
  F_prev <- c(0, Fa[-m])
  denom <- 1 - Fa
  denom_prev <- 1 - F_prev
  R <- numeric(m)
  for (k in seq_len(m)) {
    risk_after <- as.numeric(y_i > timepoints[k])
    risk_before <- if (k == 1L) 1 else as.numeric(y_i > timepoints[k - 1L])
    cens_at <- as.numeric(delta_i == 0 && y_i == timepoints[k])
    R[k] <- (risk_after + cens_at) / denom[k] - risk_before / denom_prev[k]
  }
  R
}

cir_sigma.cir <- function(Q, Fa, G_i) {
  m <- length(Fa)
  F_prev <- c(0, Fa[-m])
  sig <- (1 - G_i) * (1 / (1 - Fa) - 1 / (1 - F_prev))
  Q %*% diag(sig, m, m) %*% t(Q)
}

cir_score_matrix <- function(rho, a, fun_tx, F0_nuis, e_vec, G_mat, y, delta_c,
                             timepoints = NULL, rho.weight = rho) {
  n <- length(a)
  m <- ncol(F0_nuis)
  p <- length(rho)
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  # The estimating residual Z_a(rho) is evaluated at the candidate rho, while
  # the efficient outer weight is the plug-in weight at the nuisance estimate
  # rho.weight (Remark A4 in the paper).
  F_score <- candidate_cif.cir(rho, fun_tx, F0_nuis)
  F_weight <- candidate_cif.cir(rho.weight, fun_tx, F0_nuis)
  S_mat <- matrix(0, n, p)

  for (i in seq_len(n)) {
    F_i <- fun_tx_row_stack.cir(fun_tx, i, m)
    e_i <- e_vec[i]
    Q0_score <- cir_build_Q.cir(F_score$F0[i, ])
    Q1_score <- cir_build_Q.cir(F_score$F1[i, ])
    R0 <- cir_R.cir(F_score$F0[i, ], y[i], delta_c[i], timepoints)
    R1 <- cir_R.cir(F_score$F1[i, ], y[i], delta_c[i], timepoints)
    Sigma0_score <- cir_sigma.cir(Q0_score, F_score$F0[i, ], G_mat[i, ])
    Sigma1_score <- cir_sigma.cir(Q1_score, F_score$F1[i, ], G_mat[i, ])
    J0 <- Q0_score %*% R0
    J1 <- Q1_score %*% R1
    Z0 <- as.numeric(safe_solve.cir(Sigma0_score, J0))
    Z1 <- as.numeric(safe_solve.cir(Sigma1_score, J1))
    contrast <- (a[i] / e_i) * Z1 - ((1 - a[i]) / (1 - e_i)) * Z0

    Q0_weight <- cir_build_Q.cir(F_weight$F0[i, ])
    Q1_weight <- cir_build_Q.cir(F_weight$F1[i, ])
    Sigma0_weight <- cir_sigma.cir(Q0_weight, F_weight$F0[i, ], G_mat[i, ])
    Sigma1_weight <- cir_sigma.cir(Q1_weight, F_weight$F1[i, ], G_mat[i, ])
    SigmaX_weight <- e_i * Sigma1_weight + (1 - e_i) * Sigma0_weight
    Omega <- e_i * (1 - e_i) * Sigma1_weight %*%
      safe_solve.cir(SigmaX_weight, Sigma0_weight)
    S_mat[i, ] <- as.vector(t(F_i) %*% Omega %*% contrast)
  }
  S_mat
}

dr.estimate.cir <- function(a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
                            rho.nuis, tau.nuis, gamma.nuis, eta.nuis,
                            time.dep = FALSE, rho.start = NULL,
                            max.step = 1000, thres = 1e-6, timepoints = NULL) {
  x_sop <- as_matrix.cir(x_sop, n = length(a), name = "x_sop")
  x_ps <- as_matrix.cir(x_ps, n = length(a), name = "x_ps")
  x_cen <- as_matrix.cir(x_cen, n = length(a), name = "x_cen")
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  n <- length(a)
  m <- length(timepoints)
  p <- fun_tx_dim.cir(fun_tx)
  if (is.null(rho.start)) rho.start <- rho.nuis

  F_nuis <- estimate_cif0.cir(rho.nuis, tau.nuis, fun_tx, x_sop, timepoints)
  e_vec <- expit.cir(drop(x_ps %*% gamma.nuis))
  G_mat <- predict_censoring_gG.cir(y, x_cen, eta.nuis, time.dep, timepoints)$G

  objective <- function(rho) {
    S_mat <- tryCatch({
      cir_score_matrix(
        rho, a, fun_tx, F_nuis$F0, e_vec, G_mat, y, delta_c, timepoints,
        rho.weight = rho.nuis
      )
    }, error = function(e) NULL)
    if (is.null(S_mat) || any(!is.finite(S_mat))) {
      return(1e20)
    }
    mean_score <- colMeans(S_mat)
    sum(mean_score^2)
  }

  opt <- stats::optim(
    rho.start, objective,
    control = list(maxit = max.step, reltol = thres)
  )
  S_mat <- cir_score_matrix(
    opt$par, a, fun_tx, F_nuis$F0, e_vec, G_mat, y, delta_c, timepoints,
    rho.weight = rho.nuis
  )
  list(
    point.esti = opt$par,
    convergence = opt$convergence,
    value = opt$value,
    S_mat = S_mat,
    F0_nuis = F_nuis$F0,
    e_vec = e_vec,
    G_mat = G_mat,
    timepoints = timepoints
  )
}

DREst.cir <- function(a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
                      rho.nuis = NULL, tau.nuis = NULL, gamma.nuis = NULL,
                      eta.nuis = NULL, rho.start = NULL, tau.start = NULL,
                      eta.start = NULL, time.dep = FALSE,
                      max.step = 1000, thres = 1e-6, sandwich = TRUE,
                      timepoints = NULL) {
  x_sop <- as_matrix.cir(x_sop, n = length(a), name = "x_sop")
  x_ps <- as_matrix.cir(x_ps, n = length(a), name = "x_ps")
  x_cen <- as_matrix.cir(x_cen, n = length(a), name = "x_cen")
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  p_rho <- fun_tx_dim.cir(fun_tx)
  p_tau <- ncol(x_sop)

  mle <- NULL
  if (is.null(rho.nuis) || is.null(tau.nuis)) {
    mle <- MLEst.cir(a, fun_tx, x_sop, y, delta_c, rho.start, tau.start, max.step, thres, timepoints)
    rho.nuis <- mle$point.est[seq_len(p_rho)]
    tau.nuis <- mle$point.est[p_rho + seq_len(p_tau)]
  }
  ps <- NULL
  if (is.null(gamma.nuis)) {
    ps <- fit_ps.cir(a, x_ps)
    gamma.nuis <- ps$point.est
  }
  cen <- NULL
  if (is.null(eta.nuis)) {
    cen <- fit_censoring.cir(y, delta_c, x_cen, eta.start, max.step, thres, time.dep, timepoints)
    eta.nuis <- cen$point.est
  }

  dr <- dr.estimate.cir(
    a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
    rho.nuis, tau.nuis, gamma.nuis, eta.nuis,
    time.dep = time.dep,
    rho.start = if (is.null(rho.start)) rho.nuis else rho.start,
    max.step = max.step, thres = thres, timepoints = timepoints
  )

  if (sandwich) {
    if (is.null(mle)) {
      mle <- MLEst.cir(a, fun_tx, x_sop, y, delta_c, rho.nuis, tau.nuis, max.step, thres, timepoints)
    }
    if (is.null(ps)) ps <- fit_ps.cir(a, x_ps)
    if (is.null(cen)) cen <- fit_censoring.cir(y, delta_c, x_cen, eta.nuis, max.step, thres, time.dep, timepoints)
    var <- var.DR.sandwich.cir(
      dr.point = dr$point.esti,
      a = a, fun_tx = fun_tx, x_sop = x_sop, x_ps = x_ps, x_cen = x_cen,
      y = y, delta_c = delta_c,
      rho.nuis = rho.nuis, tau.nuis = tau.nuis,
      gamma.nuis = gamma.nuis, eta.nuis = eta.nuis,
      S_mat_dr = dr$S_mat, mle = mle, ps = ps, cen = cen,
      time.dep = time.dep, timepoints = timepoints
    )
    sd.est <- var$sd.est
    cov.mat <- var$cov.mat
  } else {
    cov.mat <- safe_solve.cir(t(dr$S_mat) %*% dr$S_mat)
    sd.est <- sqrt(diag(cov.mat))
    var <- NULL
  }

  conf.lower <- dr$point.esti + stats::qnorm(0.025) * sd.est
  conf.upper <- dr$point.esti + stats::qnorm(0.975) * sd.est
  p.temp <- stats::pnorm(dr$point.esti / sd.est, 0, 1)
  p.value <- 2 * ifelse(p.temp < 0.5, p.temp, 1 - p.temp)
  names(dr$point.esti) <- names(sd.est) <- names(conf.lower) <- names(conf.upper) <-
    names(p.value) <- paste0("rho", seq_len(p_rho) - 1L)

  list(
    dr.point = dr$point.esti,
    sd.est = sd.est,
    cov.mat = cov.mat,
    conf.lower = conf.lower,
    conf.upper = conf.upper,
    p.value = p.value,
    convergence = dr$convergence,
    S_mat = dr$S_mat,
    mle = mle,
    ps = ps,
    cen = cen,
    sandwich = var
  )
}
