################################################################################
# Analytic derivatives and sandwich variance for the CIR estimator
#
# Source this file AFTER:
#   helpers_cir.R, getProb_sequent_cir.R, MLE_point_cir.R, CallMLE_cir.R,
#   nuisance_cir_new.R, DR_point_cir.R
#
# The functions below do not use numDeriv. They leave the original functions
# unchanged and add *.analytic.cir versions.
################################################################################

cir_or_null.cir <- function(x, y) if (is.null(x)) y else x

################################################################################
# 1. Derivatives of the RR + odds-product inverse map
################################################################################

cir_rr_logprob_derivatives.cir <- function(q0, q1) {
  if (length(q0) != 1L || length(q1) != 1L ||
      !is.finite(q0) || !is.finite(q1) ||
      q0 <= 0 || q0 >= 1 || q1 <= 0 || q1 >= 1) {
    stop("q0 and q1 must be finite scalar probabilities in (0, 1).",
         call. = FALSE)
  }

  s0 <- 1 - q0
  s1 <- 1 - q1
  den <- s0 + s1

  ## Gamma(psi, b) = log(q0), where log(q1/q0) = psi and log(OP) = b.
  Gamma_psi <- -s0 / den
  Gamma_b <- s0 * s1 / den

  q0_psi <- q0 * Gamma_psi
  q1_psi <- q1 * (Gamma_psi + 1)
  q0_b <- q0 * Gamma_b
  q1_b <- q1 * Gamma_b

  d_Gamma_psi <- function(dq0, dq1) {
    ds0 <- -dq0
    ds1 <- -dq1
    dden <- ds0 + ds1
    -(ds0 * den - s0 * dden) / den^2
  }

  d_Gamma_b <- function(dq0, dq1) {
    ds0 <- -dq0
    ds1 <- -dq1
    dden <- ds0 + ds1
    ((ds0 * s1 + s0 * ds1) * den - s0 * s1 * dden) / den^2
  }

  Gamma_psipsi <- d_Gamma_psi(q0_psi, q1_psi)
  Gamma_psib <- d_Gamma_psi(q0_b, q1_b)
  Gamma_bpsi <- d_Gamma_b(q0_psi, q1_psi)
  Gamma_bb <- d_Gamma_b(q0_b, q1_b)

  ## The two mixed derivatives are algebraically identical. Averaging removes
  ## negligible floating-point asymmetry.
  Gamma_psib <- 0.5 * (Gamma_psib + Gamma_bpsi)

  list(
    Gamma_psi = Gamma_psi,
    Gamma_b = Gamma_b,
    Gamma_psipsi = Gamma_psipsi,
    Gamma_psib = Gamma_psib,
    Gamma_bb = Gamma_bb
  )
}

################################################################################
# 2. Fully analytic CIR likelihood score and observed Hessian
################################################################################

cir_mle_components_analytic.cir <- function(
    pars, a, fun_tx, x_sop, y, delta_c, timepoints,
    need_hessian = FALSE, probability_floor = 1e-300) {

  n <- length(a)
  m <- length(timepoints)
  x_sop <- as_matrix.cir(x_sop, n = n, name = "x_sop")
  p_rho <- fun_tx_dim.cir(fun_tx)
  p_tau <- ncol(x_sop)
  p_all <- p_rho + p_tau

  rho <- pars[seq_len(p_rho)]
  tau <- pars[p_rho + seq_len(p_tau)]
  p0p1 <- getProb.sequent.cir(rho, tau, fun_tx, x_sop, timepoints)

  loglik_i <- numeric(n)
  score_i <- matrix(0, nrow = n, ncol = p_all)
  hessian_sum <- if (need_hessian) matrix(0, p_all, p_all) else NULL

  rho_idx <- seq_len(p_rho)
  tau_idx <- p_rho + seq_len(p_tau)

  for (i in seq_len(n)) {
    arm01 <- as.integer(a[i])
    arm_col <- arm01 + 1L
    xi <- as.numeric(x_sop[i, , drop = TRUE])

    q_arm <- numeric(m)
    dlogq <- matrix(0, nrow = m, ncol = p_all)
    Hlogq <- if (need_hessian) array(0, dim = c(p_all, p_all, m)) else NULL

    for (j in seq_len(m)) {
      q0 <- p0p1[[j]][i, 1L]
      q1 <- p0p1[[j]][i, 2L]
      dd <- cir_rr_logprob_derivatives.cir(q0, q1)

      dpsi <- as.numeric(
        fun_tx_matrix.cir(fun_tx, j, n = n)[i, , drop = TRUE] -
          fun_tx_matrix.cir(fun_tx, j + 1L, n = n)[i, , drop = TRUE]
      )

      q_arm[j] <- if (arm01 == 0L) q0 else q1
      dlogq[j, rho_idx] <- (dd$Gamma_psi + arm01) * dpsi
      dlogq[j, tau_idx] <- dd$Gamma_b * xi

      if (need_hessian) {
        Hj <- matrix(0, p_all, p_all)
        Hj[rho_idx, rho_idx] <- dd$Gamma_psipsi * tcrossprod(dpsi)
        cross_block <- dd$Gamma_psib * outer(dpsi, xi)
        Hj[rho_idx, tau_idx] <- cross_block
        Hj[tau_idx, rho_idx] <- t(cross_block)
        Hj[tau_idx, tau_idx] <- dd$Gamma_bb * tcrossprod(xi)
        Hlogq[, , j] <- Hj
      }
    }

    ## F_A(t_k) = product_{j=k}^m q_{A,j}.
    F_vec <- numeric(m)
    dF <- matrix(0, nrow = m, ncol = p_all)
    HF <- if (need_hessian) array(0, dim = c(p_all, p_all, m)) else NULL

    acc_log <- 0
    acc_score <- numeric(p_all)
    acc_hessian <- if (need_hessian) matrix(0, p_all, p_all) else NULL

    for (j in rev(seq_len(m))) {
      acc_log <- acc_log + log(q_arm[j])
      acc_score <- acc_score + dlogq[j, ]
      if (need_hessian) acc_hessian <- acc_hessian + Hlogq[, , j]

      F_vec[j] <- exp(acc_log)
      dF[j, ] <- F_vec[j] * acc_score
      if (need_hessian) {
        HF[, , j] <- F_vec[j] * (acc_hessian + tcrossprod(acc_score))
      }
    }

    k <- match(y[i], timepoints)
    if (is.na(k)) stop("An observed y value is absent from timepoints.", call. = FALSE)

    if (delta_c[i] == 1L) {
      F_prev <- if (k == 1L) 0 else F_vec[k - 1L]
      dF_prev <- if (k == 1L) numeric(p_all) else dF[k - 1L, ]
      HF_prev <- if (!need_hessian || k == 1L) {
        matrix(0, p_all, p_all)
      } else {
        HF[, , k - 1L]
      }

      mass <- F_vec[k] - F_prev
      if (!is.finite(mass) || mass <= 0) {
        stop("The coherent CIR likelihood produced a non-positive event mass.",
             call. = FALSE)
      }
      mass_eval <- max(mass, probability_floor)
      dmass <- dF[k, ] - dF_prev

      loglik_i[i] <- log(mass_eval)
      score_i[i, ] <- dmass / mass_eval
      if (need_hessian) {
        Hmass <- HF[, , k] - HF_prev
        hessian_sum <- hessian_sum +
          Hmass / mass_eval - tcrossprod(dmass) / mass_eval^2
      }
    } else {
      surv <- 1 - F_vec[k]
      if (!is.finite(surv) || surv <= 0) {
        stop("The coherent CIR likelihood produced a non-positive survival probability.",
             call. = FALSE)
      }
      surv_eval <- max(surv, probability_floor)

      loglik_i[i] <- log(surv_eval)
      score_i[i, ] <- -dF[k, ] / surv_eval
      if (need_hessian) {
        hessian_sum <- hessian_sum -
          HF[, , k] / surv_eval - tcrossprod(dF[k, ]) / surv_eval^2
      }
    }
  }

  list(
    value = -sum(loglik_i),
    gradient = -colSums(score_i),
    hessian = if (need_hessian) -hessian_sum else NULL,
    loglik_i = loglik_i,
    S_mat = score_i
  )
}

max.likelihood.cir.analytic <- function(
    a, fun_tx, x_sop, y, delta_c, rho.start = NULL, tau.start = NULL,
    max.step = 1000, thres = 1e-8, timepoints = NULL,
    vcov.type = c("model", "sandwich", "opg")) {

  vcov.type <- match.arg(vcov.type)
  x_sop <- as_matrix.cir(x_sop, n = length(a), name = "x_sop")
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  p_rho <- fun_tx_dim.cir(fun_tx)
  p_tau <- ncol(x_sop)
  if (is.null(rho.start)) rho.start <- rep(0, p_rho)
  if (is.null(tau.start)) tau.start <- rep(0, p_tau)
  start <- c(rho.start, tau.start)

  cache <- new.env(parent = emptyenv())
  cache$par <- NULL
  cache$value <- NULL
  evaluate <- function(pars) {
    same <- !is.null(cache$par) && length(cache$par) == length(pars) &&
      isTRUE(all(pars == cache$par))
    if (!same) {
      cache$par <- pars
      cache$value <- cir_mle_components_analytic.cir(
        pars, a, fun_tx, x_sop, y, delta_c, timepoints,
        need_hessian = FALSE
      )
    }
    cache$value
  }

  opt <- stats::optim(
    par = start,
    fn = function(pars) evaluate(pars)$value,
    gr = function(pars) evaluate(pars)$gradient,
    method = "BFGS",
    control = list(maxit = max.step, reltol = thres)
  )

  final <- cir_mle_components_analytic.cir(
    opt$par, a, fun_tx, x_sop, y, delta_c, timepoints,
    need_hessian = TRUE
  )

  bread <- safe_solve.cir(final$hessian)
  meat <- crossprod(final$S_mat)
  cov.model <- bread
  cov.sandwich <- bread %*% meat %*% t(bread)
  cov.opg <- safe_solve.cir(meat)
  cov.use <- switch(
    vcov.type,
    model = cov.model,
    sandwich = cov.sandwich,
    opg = cov.opg
  )

  list(
    par = opt$par,
    convergence = opt$convergence,
    value = final$value,
    cov = cov.use,
    bread = bread,
    cov.model = cov.model,
    cov.sandwich = cov.sandwich,
    cov.opg = cov.opg,
    S_mat = final$S_mat,
    hessian = final$hessian,
    gradient = final$gradient,
    counts = opt$counts,
    message = opt$message,
    timepoints = timepoints,
    vcov.type = vcov.type
  )
}

MLEst.cir.analytic <- function(
    a, fun_tx, x_sop, y, delta_c, rho.start = NULL, tau.start = NULL,
    max.step = 1000, thres = 1e-8, timepoints = NULL,
    vcov.type = c("model", "sandwich", "opg"), conf.level = 0.95) {

  vcov.type <- match.arg(vcov.type)
  x_sop <- as_matrix.cir(x_sop, n = length(a), name = "x_sop")
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  check_fun_tx.cir(fun_tx, n = length(a), m = length(timepoints))

  p_rho <- fun_tx_dim.cir(fun_tx)
  p_tau <- ncol(x_sop)
  fit <- max.likelihood.cir.analytic(
    a, fun_tx, x_sop, y, delta_c,
    rho.start = rho.start, tau.start = tau.start,
    max.step = max.step, thres = thres, timepoints = timepoints,
    vcov.type = vcov.type
  )

  point.est <- fit$par
  se.est <- sqrt(pmax(diag(fit$cov), 0))
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  nm <- c(paste0("rho", seq_len(p_rho) - 1L),
          paste0("tau", seq_len(p_tau) - 1L))

  names(point.est) <- names(se.est) <- nm
  rownames(fit$cov) <- colnames(fit$cov) <- nm
  rownames(fit$bread) <- colnames(fit$bread) <- nm

  list(
    point.est = point.est,
    se.est = se.est,
    cov = fit$cov,
    bread = fit$bread,
    cov.model = fit$cov.model,
    cov.sandwich = fit$cov.sandwich,
    cov.opg = fit$cov.opg,
    S_mat = fit$S_mat,
    hessian = fit$hessian,
    conf.lower = point.est - z * se.est,
    conf.upper = point.est + z * se.est,
    p.value = 2 * stats::pnorm(-abs(point.est / se.est)),
    convergence = fit$convergence,
    gradient = fit$gradient,
    counts = fit$counts,
    message = fit$message,
    timepoints = fit$timepoints,
    vcov.type = vcov.type
  )
}

################################################################################
# 3. Analytic censoring score, Hessian, and derivative of G(t_k; X, eta)
################################################################################

cir_censor_design_row.cir <- function(i, j, x_cen, m, time.dep) {
  if (!time.dep) return(as.numeric(x_cen[i, , drop = TRUE]))

  p_cov <- ncol(x_cen) - 1L
  z <- numeric((m - 1L) + p_cov)
  z[j] <- 1
  if (p_cov > 0L) z[(m):(m - 1L + p_cov)] <- x_cen[i, -1L, drop = TRUE]
  z
}

censoring_components_analytic.cir <- function(
    eta, y, delta_c, x_cen, timepoints, time.dep = FALSE,
    need_hessian = TRUE, probability_floor = 1e-12) {

  n <- length(y)
  m <- length(timepoints)
  x_cen <- as_matrix.cir(x_cen, n = n, name = "x_cen")
  p_eta <- length(eta)

  if (m <= 1L) {
    return(list(
      value = 0,
      gradient = rep(0, p_eta),
      hessian = matrix(0, p_eta, p_eta),
      S_mat = matrix(0, n, p_eta),
      loglik_i = rep(0, n)
    ))
  }

  lambda <- censoring_hazard.cir(x_cen, eta, m, time.dep)
  S_mat <- matrix(0, nrow = n, ncol = p_eta)
  Hneg <- matrix(0, p_eta, p_eta)
  loglik_i <- numeric(n)

  for (i in seq_len(n)) {
    k_obs <- match(y[i], timepoints)
    if (is.na(k_obs)) stop("An observed y value is absent from timepoints.", call. = FALSE)

    for (j in seq_len(m - 1L)) {
      cens_event <- as.numeric(delta_c[i] == 0L && k_obs == j)
      observed_hazard <- as.numeric(k_obs > j || cens_event == 1)
      if (observed_hazard == 0) next

      lam <- min(max(lambda[i, j], probability_floor), 1 - probability_floor)
      z <- cir_censor_design_row.cir(i, j, x_cen, m, time.dep)

      loglik_i[i] <- loglik_i[i] +
        cens_event * log(lam) + (observed_hazard - cens_event) * log1p(-lam)
      S_mat[i, ] <- S_mat[i, ] + z * (cens_event - observed_hazard * lam)
      if (need_hessian) {
        Hneg <- Hneg + observed_hazard * lam * (1 - lam) * tcrossprod(z)
      }
    }
  }

  list(
    value = -sum(loglik_i),
    gradient = -colSums(S_mat),
    hessian = if (need_hessian) Hneg else NULL,
    S_mat = S_mat,
    loglik_i = loglik_i
  )
}

fit_censoring.cir.analytic <- function(
    y, delta_c, x_cen, eta.start = NULL, max.step = 1000, thres = 1e-8,
    time.dep = FALSE, timepoints = NULL) {

  x_cen <- as_matrix.cir(x_cen, n = length(y), name = "x_cen")
  timepoints <- resolve_timepoints.cir(y, timepoints = timepoints)
  m <- length(timepoints)

  if (is.null(eta.start)) {
    eta.start <- if (time.dep) {
      rep(0, (m - 1L) + ncol(x_cen) - 1L)
    } else {
      rep(0, ncol(x_cen))
    }
  }

  cache <- new.env(parent = emptyenv())
  cache$par <- NULL
  cache$value <- NULL
  evaluate <- function(pars) {
    same <- !is.null(cache$par) && length(cache$par) == length(pars) &&
      isTRUE(all(pars == cache$par))
    if (!same) {
      cache$par <- pars
      cache$value <- censoring_components_analytic.cir(
        pars, y, delta_c, x_cen, timepoints,
        time.dep = time.dep, need_hessian = FALSE
      )
    }
    cache$value
  }

  opt <- stats::optim(
    par = eta.start,
    fn = function(pars) evaluate(pars)$value,
    gr = function(pars) evaluate(pars)$gradient,
    method = "BFGS",
    control = list(maxit = max.step, reltol = thres)
  )

  final <- censoring_components_analytic.cir(
    opt$par, y, delta_c, x_cen, timepoints,
    time.dep = time.dep, need_hessian = TRUE
  )
  bread <- safe_solve.cir(final$hessian)

  list(
    point.est = opt$par,
    cov = bread,
    bread = bread,
    S_mat = final$S_mat,
    hessian = final$hessian,
    convergence = opt$convergence,
    value = final$value,
    gradient = final$gradient,
    counts = opt$counts,
    message = opt$message,
    timepoints = timepoints,
    time.dep = time.dep
  )
}

predict_censoring_G_derivative.cir <- function(
    y, x_cen, eta, time.dep = FALSE, timepoints = NULL) {

  x_cen <- as_matrix.cir(x_cen, n = length(y), name = "x_cen")
  timepoints <- resolve_timepoints.cir(y, timepoints = timepoints)
  n <- nrow(x_cen)
  m <- length(timepoints)
  p_eta <- length(eta)

  if (m <= 1L) {
    return(list(
      G = matrix(0, n, 1L),
      dG = array(0, dim = c(n, 1L, p_eta)),
      lambda = matrix(numeric(0), nrow = n, ncol = 0)
    ))
  }

  lambda <- censoring_hazard.cir(x_cen, eta, m, time.dep)
  G <- matrix(0, nrow = n, ncol = m)
  dG <- array(0, dim = c(n, m, p_eta))

  for (i in seq_len(n)) {
    surv <- 1
    cumulative_direction <- numeric(p_eta)
    for (j in seq_len(m - 1L)) {
      z <- cir_censor_design_row.cir(i, j, x_cen, m, time.dep)
      lam <- lambda[i, j]
      surv <- surv * (1 - lam)
      cumulative_direction <- cumulative_direction + lam * z
      G[i, j + 1L] <- 1 - surv
      dG[i, j + 1L, ] <- surv * cumulative_direction
    }
  }

  list(G = G, dG = dG, lambda = lambda)
}

################################################################################
# 4. First derivative of nuisance F0(t_k; X, rho_nuis, tau_nuis)
################################################################################

cir_F0_nuisance_derivative.cir <- function(
    rho.nuis, tau.nuis, fun_tx, x_sop, timepoints) {

  x_sop <- as_matrix.cir(x_sop, name = "x_sop")
  n <- nrow(x_sop)
  m <- length(timepoints)
  p_rho <- length(rho.nuis)
  p_tau <- length(tau.nuis)
  p_out <- p_rho + p_tau

  p0p1 <- getProb.sequent.cir(
    rho.nuis, tau.nuis, fun_tx, x_sop, timepoints
  )
  cif <- cif_from_sequential.cir(p0p1)
  F0 <- do.call(cbind, lapply(cif, function(z) z[, 1L]))
  F1 <- do.call(cbind, lapply(cif, function(z) z[, 2L]))
  dF0 <- array(0, dim = c(n, m, p_out))

  for (i in seq_len(n)) {
    accum <- numeric(p_out)
    xi <- as.numeric(x_sop[i, , drop = TRUE])

    for (j in rev(seq_len(m))) {
      q0 <- p0p1[[j]][i, 1L]
      q1 <- p0p1[[j]][i, 2L]
      dd <- cir_rr_logprob_derivatives.cir(q0, q1)
      dpsi <- as.numeric(
        fun_tx_matrix.cir(fun_tx, j, n = n)[i, , drop = TRUE] -
          fun_tx_matrix.cir(fun_tx, j + 1L, n = n)[i, , drop = TRUE]
      )
      vj <- c(dd$Gamma_psi * dpsi, dd$Gamma_b * xi)
      accum <- accum + vj
      dF0[i, j, ] <- F0[i, j] * accum
    }
  }

  list(F0 = F0, F1 = F1, dF0 = dF0, p0p1 = p0p1)
}

################################################################################
# 5. Directional derivatives of Q, R, Sigma, Z, and Omega
################################################################################

cir_shift_matrix.cir <- function(m) {
  Gamma <- matrix(0, m, m)
  if (m > 1L) Gamma[cbind(seq_len(m - 1L), 2:m)] <- 1
  Gamma
}

cir_assert_cdf.cir <- function(Fa, label = "F") {
  if (any(!is.finite(Fa)) || any(Fa <= 0) || any(Fa >= 1) ||
      any(diff(c(0, Fa)) <= 0)) {
    stop(label, " is not a strictly increasing CDF in (0, 1). ",
         "The DR search has left the admissible parameter region.",
         call. = FALSE)
  }
  invisible(TRUE)
}

cir_Q_direction.cir <- function(Fa, dFa = NULL) {
  Fa <- as.numeric(Fa)
  m <- length(Fa)
  if (is.null(dFa)) dFa <- numeric(m)
  dFa <- as.numeric(dFa)

  F_prev <- c(0, Fa[-m])
  F_next <- c(Fa[-1L], 1)
  dF_prev <- c(0, dFa[-m])
  dF_next <- c(dFa[-1L], 0)

  den1 <- Fa - F_prev
  den2 <- F_next - Fa
  if (any(den1 <= 0) || any(den2 <= 0)) {
    stop("Q(F) has a non-positive CDF increment.", call. = FALSE)
  }

  N1 <- Fa * (1 - F_prev)
  dN1 <- dFa * (1 - F_prev) - Fa * dF_prev
  dden1 <- dFa - dF_prev
  q1 <- -N1 / den1
  dq1 <- -(dN1 * den1 - N1 * dden1) / den1^2

  N2 <- Fa * (1 - F_next)
  dN2 <- dFa * (1 - F_next) - Fa * dF_next
  dden2 <- dF_next - dFa
  q2 <- N2 / den2
  dq2 <- (dN2 * den2 - N2 * dden2) / den2^2
  q2[m] <- 0
  dq2[m] <- 0

  Gamma <- cir_shift_matrix.cir(m)
  Q <- diag(q1, nrow = m, ncol = m) +
    diag(q2, nrow = m, ncol = m) %*% Gamma
  dQ <- diag(dq1, nrow = m, ncol = m) +
    diag(dq2, nrow = m, ncol = m) %*% Gamma

  list(Q = Q, dQ = dQ)
}

cir_R_direction.cir <- function(
    Fa, dFa, y_i, delta_i, timepoints) {

  Fa <- as.numeric(Fa)
  dFa <- as.numeric(dFa)
  m <- length(Fa)
  F_prev <- c(0, Fa[-m])
  dF_prev <- c(0, dFa[-m])

  risk_after <- as.numeric(y_i > timepoints)
  risk_before <- c(1, as.numeric(y_i > timepoints[-m]))
  cens_at <- as.numeric(delta_i == 0L & y_i == timepoints)
  first_coeff <- risk_after + cens_at

  R <- first_coeff / (1 - Fa) - risk_before / (1 - F_prev)
  dR <- first_coeff * dFa / (1 - Fa)^2 -
    risk_before * dF_prev / (1 - F_prev)^2

  list(R = R, dR = dR)
}

cir_sigma_direction.cir <- function(
    Q, dQ, Fa, dFa, G_i, dG_i) {

  Fa <- as.numeric(Fa)
  dFa <- as.numeric(dFa)
  G_i <- as.numeric(G_i)
  dG_i <- as.numeric(dG_i)
  m <- length(Fa)

  F_prev <- c(0, Fa[-m])
  dF_prev <- c(0, dFa[-m])
  diff_inv <- 1 / (1 - Fa) - 1 / (1 - F_prev)
  ddiff_inv <- dFa / (1 - Fa)^2 - dF_prev / (1 - F_prev)^2

  sigma_diag <- (1 - G_i) * diff_inv
  dsigma_diag <- -dG_i * diff_inv + (1 - G_i) * ddiff_inv
  D <- diag(sigma_diag, nrow = m, ncol = m)
  dD <- diag(dsigma_diag, nrow = m, ncol = m)

  Sigma <- Q %*% D %*% t(Q)
  dSigma <- dQ %*% D %*% t(Q) +
    Q %*% dD %*% t(Q) +
    Q %*% D %*% t(dQ)

  list(Sigma = Sigma, dSigma = dSigma)
}

cir_arm_components_direction.cir <- function(
    Fa, G_i, y_i, delta_i, timepoints,
    dFa = NULL, dG_i = NULL) {

  m <- length(Fa)
  if (is.null(dFa)) dFa <- numeric(m)
  if (is.null(dG_i)) dG_i <- numeric(m)

  q <- cir_Q_direction.cir(Fa, dFa)
  r <- cir_R_direction.cir(Fa, dFa, y_i, delta_i, timepoints)
  s <- cir_sigma_direction.cir(q$Q, q$dQ, Fa, dFa, G_i, dG_i)

  J <- as.numeric(q$Q %*% r$R)
  dJ <- as.numeric(q$dQ %*% r$R + q$Q %*% r$dR)
  Z <- as.numeric(safe_solve.cir(s$Sigma, J))
  dZ <- as.numeric(safe_solve.cir(s$Sigma, dJ - s$dSigma %*% Z))

  list(
    Q = q$Q,
    dQ = q$dQ,
    R = r$R,
    dR = r$dR,
    Sigma = s$Sigma,
    dSigma = s$dSigma,
    J = J,
    dJ = dJ,
    Z = Z,
    dZ = dZ
  )
}

cir_omega_direction.cir <- function(
    e, Sigma1, Sigma0, de = 0,
    dSigma1 = NULL, dSigma0 = NULL) {

  m <- nrow(Sigma0)
  if (is.null(dSigma1)) dSigma1 <- matrix(0, m, m)
  if (is.null(dSigma0)) dSigma0 <- matrix(0, m, m)

  SigmaX <- e * Sigma1 + (1 - e) * Sigma0
  K <- safe_solve.cir(SigmaX, Sigma0)
  coef <- e * (1 - e)
  Omega <- coef * Sigma1 %*% K

  dSigmaX <- de * (Sigma1 - Sigma0) +
    e * dSigma1 + (1 - e) * dSigma0
  dK <- safe_solve.cir(SigmaX, dSigma0 - dSigmaX %*% K)
  dcoef <- (1 - 2 * e) * de
  dOmega <- dcoef * Sigma1 %*% K +
    coef * (dSigma1 %*% K + Sigma1 %*% dK)

  list(Omega = Omega, dOmega = dOmega, SigmaX = SigmaX)
}

################################################################################
# 6. Analytic Jacobians of the full DR estimating equation
################################################################################

cir_score_jacobians_analytic.cir <- function(
    dr.point, a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
    rho.nuis, tau.nuis, gamma.nuis, eta.nuis,
    time.dep = FALSE, timepoints = NULL, check_cdf = TRUE) {

  n <- length(a)
  x_sop <- as_matrix.cir(x_sop, n = n, name = "x_sop")
  x_ps <- as_matrix.cir(x_ps, n = n, name = "x_ps")
  x_cen <- as_matrix.cir(x_cen, n = n, name = "x_cen")
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  m <- length(timepoints)

  p_target <- length(dr.point)
  p_rho_nuis <- length(rho.nuis)
  p_tau <- length(tau.nuis)
  p_gamma <- length(gamma.nuis)
  p_eta <- length(eta.nuis)
  p_outcome <- p_rho_nuis + p_tau
  p_nuisance <- p_outcome + p_gamma + p_eta

  if (p_target != p_rho_nuis) {
    stop("dr.point and rho.nuis must have the same dimension.", call. = FALSE)
  }

  f0fit <- cir_F0_nuisance_derivative.cir(
    rho.nuis, tau.nuis, fun_tx, x_sop, timepoints
  )
  e_vec <- expit.cir(drop(x_ps %*% gamma.nuis))
  cen_deriv <- predict_censoring_G_derivative.cir(
    y, x_cen, eta.nuis, time.dep = time.dep, timepoints = timepoints
  )
  G_mat <- cen_deriv$G

  S_mat <- matrix(0, nrow = n, ncol = p_target)
  dU.by.rho <- matrix(0, nrow = p_target, ncol = p_target)
  dU.by.nuis <- matrix(0, nrow = p_target, ncol = p_nuisance)

  out_cols <- seq_len(p_outcome)
  gamma_cols <- p_outcome + seq_len(p_gamma)
  eta_cols <- p_outcome + p_gamma + seq_len(p_eta)

  for (i in seq_len(n)) {
    F_design <- fun_tx_row_stack.cir(fun_tx, i, m)
    F0 <- as.numeric(f0fit$F0[i, ])
    eta_score <- as.numeric(F_design %*% dr.point)
    eta_weight <- as.numeric(F_design %*% rho.nuis)
    F1_score <- F0 * exp(eta_score)
    F1_weight <- F0 * exp(eta_weight)
    G_i <- as.numeric(G_mat[i, ])
    e_i <- e_vec[i]

    if (check_cdf) {
      cir_assert_cdf.cir(F0, paste0("F0 for observation ", i))
      cir_assert_cdf.cir(F1_score, paste0("candidate F1 for observation ", i))
      cir_assert_cdf.cir(F1_weight, paste0("weight F1 for observation ", i))
      if (!is.finite(e_i) || e_i <= 0 || e_i >= 1) {
        stop("Propensity score violates positivity for observation ", i, ".",
             call. = FALSE)
      }
    }

    zero_m <- numeric(m)
    arm0 <- cir_arm_components_direction.cir(
      F0, G_i, y[i], delta_c[i], timepoints, zero_m, zero_m
    )
    arm1_score <- cir_arm_components_direction.cir(
      F1_score, G_i, y[i], delta_c[i], timepoints, zero_m, zero_m
    )
    arm1_weight <- cir_arm_components_direction.cir(
      F1_weight, G_i, y[i], delta_c[i], timepoints, zero_m, zero_m
    )

    omega0 <- cir_omega_direction.cir(
      e_i, arm1_weight$Sigma, arm0$Sigma
    )
    Omega <- omega0$Omega
    contrast <- (a[i] / e_i) * arm1_score$Z -
      ((1 - a[i]) / (1 - e_i)) * arm0$Z
    S_mat[i, ] <- as.numeric(crossprod(F_design, Omega %*% contrast))

    ## Derivative with respect to the target rho. The plug-in efficient weight
    ## is held fixed at rho.nuis, exactly as in Remark A4.
    for (ell in seq_len(p_target)) {
      dF1 <- F1_score * F_design[, ell]
      d_arm1 <- cir_arm_components_direction.cir(
        F1_score, G_i, y[i], delta_c[i], timepoints,
        dFa = dF1, dG_i = zero_m
      )
      dcontrast <- (a[i] / e_i) * d_arm1$dZ
      dUi <- as.numeric(crossprod(F_design, Omega %*% dcontrast))
      dU.by.rho[, ell] <- dU.by.rho[, ell] + dUi
    }

    ## Derivative with respect to outcome nuisance (rho.nuis, tau.nuis).
    for (ell in seq_len(p_outcome)) {
      dF0 <- as.numeric(f0fit$dF0[i, , ell])
      dF1_score <- exp(eta_score) * dF0
      if (ell <= p_rho_nuis) {
        dF1_weight <- exp(eta_weight) *
          (dF0 + F0 * F_design[, ell])
      } else {
        dF1_weight <- exp(eta_weight) * dF0
      }

      d_arm0 <- cir_arm_components_direction.cir(
        F0, G_i, y[i], delta_c[i], timepoints,
        dFa = dF0, dG_i = zero_m
      )
      d_arm1_score <- cir_arm_components_direction.cir(
        F1_score, G_i, y[i], delta_c[i], timepoints,
        dFa = dF1_score, dG_i = zero_m
      )
      d_arm1_weight <- cir_arm_components_direction.cir(
        F1_weight, G_i, y[i], delta_c[i], timepoints,
        dFa = dF1_weight, dG_i = zero_m
      )
      d_omega <- cir_omega_direction.cir(
        e_i, arm1_weight$Sigma, arm0$Sigma,
        de = 0,
        dSigma1 = d_arm1_weight$dSigma,
        dSigma0 = d_arm0$dSigma
      )$dOmega

      dcontrast <- (a[i] / e_i) * d_arm1_score$dZ -
        ((1 - a[i]) / (1 - e_i)) * d_arm0$dZ
      dUi <- as.numeric(crossprod(
        F_design,
        d_omega %*% contrast + Omega %*% dcontrast
      ))
      dU.by.nuis[, out_cols[ell]] <-
        dU.by.nuis[, out_cols[ell]] + dUi
    }

    ## Derivative with respect to propensity-score nuisance gamma.
    for (ell in seq_len(p_gamma)) {
      de <- e_i * (1 - e_i) * x_ps[i, ell]
      d_omega <- cir_omega_direction.cir(
        e_i, arm1_weight$Sigma, arm0$Sigma, de = de
      )$dOmega
      dcontrast <- -a[i] * de / e_i^2 * arm1_score$Z -
        (1 - a[i]) * de / (1 - e_i)^2 * arm0$Z
      dUi <- as.numeric(crossprod(
        F_design,
        d_omega %*% contrast + Omega %*% dcontrast
      ))
      dU.by.nuis[, gamma_cols[ell]] <-
        dU.by.nuis[, gamma_cols[ell]] + dUi
    }

    ## Derivative with respect to censoring nuisance eta.
    for (ell in seq_len(p_eta)) {
      dG <- as.numeric(cen_deriv$dG[i, , ell])
      d_arm0 <- cir_arm_components_direction.cir(
        F0, G_i, y[i], delta_c[i], timepoints,
        dFa = zero_m, dG_i = dG
      )
      d_arm1_score <- cir_arm_components_direction.cir(
        F1_score, G_i, y[i], delta_c[i], timepoints,
        dFa = zero_m, dG_i = dG
      )
      d_arm1_weight <- cir_arm_components_direction.cir(
        F1_weight, G_i, y[i], delta_c[i], timepoints,
        dFa = zero_m, dG_i = dG
      )
      d_omega <- cir_omega_direction.cir(
        e_i, arm1_weight$Sigma, arm0$Sigma,
        de = 0,
        dSigma1 = d_arm1_weight$dSigma,
        dSigma0 = d_arm0$dSigma
      )$dOmega
      dcontrast <- (a[i] / e_i) * d_arm1_score$dZ -
        ((1 - a[i]) / (1 - e_i)) * d_arm0$dZ
      dUi <- as.numeric(crossprod(
        F_design,
        d_omega %*% contrast + Omega %*% dcontrast
      ))
      dU.by.nuis[, eta_cols[ell]] <-
        dU.by.nuis[, eta_cols[ell]] + dUi
    }
  }

  colnames(dU.by.nuis) <- c(
    paste0("rho.nuis", seq_len(p_rho_nuis) - 1L),
    paste0("tau.nuis", seq_len(p_tau) - 1L),
    paste0("gamma", seq_len(p_gamma) - 1L),
    paste0("eta", seq_len(p_eta) - 1L)
  )

  list(
    S_mat = S_mat,
    dU.by.rho = dU.by.rho,
    dU.by.nuis = dU.by.nuis,
    F0_nuis = f0fit$F0,
    e_vec = e_vec,
    G_mat = G_mat
  )
}

################################################################################
# 7. Analytic full sandwich variance
################################################################################

cir_fit_bread.cir <- function(fit, label) {
  if (is.null(fit)) stop(label, " fit is required for its nuisance correction.", call. = FALSE)
  bread <- cir_or_null.cir(fit$bread, fit$cov)
  if (is.null(bread)) stop(label, " fit has neither bread nor cov.", call. = FALSE)
  bread
}

var.DR.sandwich.analytic.cir <- function(
    dr.point, a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
    rho.nuis, tau.nuis, gamma.nuis, eta.nuis,
    mle = NULL, ps = NULL, cen = NULL,
    outcome.estimated = TRUE, ps.estimated = TRUE, cen.estimated = TRUE,
    time.dep = FALSE, timepoints = NULL,
    center.meat = FALSE, hc1 = FALSE, check_cdf = TRUE,
    condition.warn = 1e10) {

  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  n <- length(a)
  p_rho <- length(rho.nuis)
  p_tau <- length(tau.nuis)
  p_gamma <- length(gamma.nuis)
  p_eta <- length(eta.nuis)

  jac <- cir_score_jacobians_analytic.cir(
    dr.point = dr.point,
    a = a, fun_tx = fun_tx, x_sop = x_sop, x_ps = x_ps, x_cen = x_cen,
    y = y, delta_c = delta_c,
    rho.nuis = rho.nuis, tau.nuis = tau.nuis,
    gamma.nuis = gamma.nuis, eta.nuis = eta.nuis,
    time.dep = time.dep, timepoints = timepoints,
    check_cdf = check_cdf
  )

  i_out <- seq_len(p_rho + p_tau)
  i_ps <- max(i_out) + seq_len(p_gamma)
  i_cen <- max(i_ps) + seq_len(p_eta)

  U_tilde <- jac$S_mat
  if (outcome.estimated) {
    bread <- cir_fit_bread.cir(mle, "Outcome MLE")
    U_tilde <- U_tilde + t(
      jac$dU.by.nuis[, i_out, drop = FALSE] %*%
        bread %*% t(mle$S_mat)
    )
  }
  if (ps.estimated) {
    bread <- cir_fit_bread.cir(ps, "Propensity-score")
    U_tilde <- U_tilde + t(
      jac$dU.by.nuis[, i_ps, drop = FALSE] %*%
        bread %*% t(ps$S_mat)
    )
  }
  if (cen.estimated) {
    bread <- cir_fit_bread.cir(cen, "Censoring")
    U_tilde <- U_tilde + t(
      jac$dU.by.nuis[, i_cen, drop = FALSE] %*%
        bread %*% t(cen$S_mat)
    )
  }

  meat_scores <- if (center.meat) {
    sweep(U_tilde, 2L, colMeans(U_tilde), "-")
  } else {
    U_tilde
  }
  meat <- crossprod(meat_scores)
  if (hc1 && n > length(dr.point)) {
    meat <- meat * n / (n - length(dr.point))
  }

  condition_A <- kappa(jac$dU.by.rho)
  if (is.finite(condition.warn) && is.finite(condition_A) &&
      condition_A > condition.warn) {
    warning(
      "The target estimating-equation Jacobian is ill-conditioned: kappa = ",
      format(condition_A, digits = 4),
      ". Rare-event Wald inference may be unstable.",
      call. = FALSE
    )
  }

  Ainv <- -safe_solve.cir(jac$dU.by.rho)
  cov.mat <- Ainv %*% meat %*% t(Ainv)
  cov.mat <- 0.5 * (cov.mat + t(cov.mat))

  list(
    sd.est = sqrt(pmax(diag(cov.mat), 0)),
    cov.mat = cov.mat,
    U_tilde = U_tilde,
    dU.by.rho = jac$dU.by.rho,
    dU.by.nuis = jac$dU.by.nuis,
    S_mat = jac$S_mat,
    Ainv = Ainv,
    meat = meat,
    condition.dU.rho = condition_A,
    center.meat = center.meat,
    hc1 = hc1
  )
}


################################################################################
# 8. Checked DR point estimator: prevents optimization outside valid CDF region
################################################################################

cir_candidate_cdf_valid.cir <- function(rho, fun_tx, F0, tol = 0) {
  cand <- candidate_cif.cir(rho, fun_tx, F0)$F1
  if (any(!is.finite(cand)) || any(cand <= tol) || any(cand >= 1 - tol)) {
    return(FALSE)
  }
  if (ncol(cand) > 1L && any(cand[, -1L, drop = FALSE] -
                              cand[, -ncol(cand), drop = FALSE] <= tol)) {
    return(FALSE)
  }
  TRUE
}

dr.estimate.cir.checked <- function(
    a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
    rho.nuis, tau.nuis, gamma.nuis, eta.nuis,
    time.dep = FALSE, rho.start = NULL,
    max.step = 1000, thres = 1e-8, timepoints = NULL,
    cdf.tol = 0, score.tol = 1e-5) {

  x_sop <- as_matrix.cir(x_sop, n = length(a), name = "x_sop")
  x_ps <- as_matrix.cir(x_ps, n = length(a), name = "x_ps")
  x_cen <- as_matrix.cir(x_cen, n = length(a), name = "x_cen")
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  if (is.null(rho.start)) rho.start <- rho.nuis

  F_nuis <- estimate_cif0.cir(rho.nuis, tau.nuis, fun_tx, x_sop, timepoints)
  e_vec <- expit.cir(drop(x_ps %*% gamma.nuis))
  G_mat <- predict_censoring_gG.cir(
    y, x_cen, eta.nuis, time.dep, timepoints
  )$G

  if (!cir_candidate_cdf_valid.cir(rho.start, fun_tx, F_nuis$F0, cdf.tol)) {
    if (cir_candidate_cdf_valid.cir(rho.nuis, fun_tx, F_nuis$F0, cdf.tol)) {
      warning("rho.start is outside the valid CDF region; using rho.nuis.",
              call. = FALSE)
      rho.start <- rho.nuis
    } else {
      stop("Neither rho.start nor rho.nuis gives a valid candidate CDF.",
           call. = FALSE)
    }
  }

  objective <- function(rho) {
    if (!cir_candidate_cdf_valid.cir(rho, fun_tx, F_nuis$F0, cdf.tol)) {
      return(1e50)
    }
    S_mat <- tryCatch(
      cir_score_matrix(
        rho, a, fun_tx, F_nuis$F0, e_vec, G_mat,
        y, delta_c, timepoints, rho.weight = rho.nuis
      ),
      error = function(e) NULL
    )
    if (is.null(S_mat) || any(!is.finite(S_mat))) return(1e50)
    mean_score <- colMeans(S_mat)
    sum(mean_score^2)
  }

  opt <- stats::optim(
    par = rho.start,
    fn = objective,
    method = "Nelder-Mead",
    control = list(maxit = max.step, reltol = thres)
  )

  if (!cir_candidate_cdf_valid.cir(opt$par, fun_tx, F_nuis$F0, cdf.tol)) {
    stop("The DR optimizer ended outside the valid CDF region.", call. = FALSE)
  }
  S_mat <- cir_score_matrix(
    opt$par, a, fun_tx, F_nuis$F0, e_vec, G_mat,
    y, delta_c, timepoints, rho.weight = rho.nuis
  )
  max_score <- max(abs(colMeans(S_mat)))
  convergence <- if (identical(opt$convergence, 0L) && max_score <= score.tol) {
    0L
  } else {
    1L
  }

  list(
    point.esti = opt$par,
    convergence = convergence,
    optim.convergence = opt$convergence,
    value = opt$value,
    S_mat = S_mat,
    max.abs.mean.score = max_score,
    F0_nuis = F_nuis$F0,
    e_vec = e_vec,
    G_mat = G_mat,
    timepoints = timepoints,
    counts = opt$counts,
    message = opt$message
  )
}

################################################################################
# 9. End-to-end wrapper using analytic nuisance derivatives and sandwich
################################################################################

DREst.cir.analytic <- function(
    a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
    rho.nuis = NULL, tau.nuis = NULL, gamma.nuis = NULL, eta.nuis = NULL,
    mle.fit = NULL, ps.fit = NULL, cen.fit = NULL,
    rho.start = NULL, tau.start = NULL, eta.start = NULL,
    time.dep = FALSE, max.step = 1000, thres = 1e-8,
    sandwich = TRUE, timepoints = NULL, mle.vcov.type = "model",
    dr.method = c("checked", "original"), supplied.nuisance.fixed = FALSE,
    center.meat = FALSE, hc1 = FALSE, check_cdf = TRUE,
    condition.warn = 1e10, cdf.tol = 0, score.tol = 1e-5,
    conf.level = 0.95) {

  if (!isTRUE(sandwich)) {
    stop(
      "DREst.cir.analytic implements the full A^{-1} B A^{-T} sandwich. ",
      "Keep sandwich = TRUE; fixed nuisance parameters are handled through ",
      "supplied.nuisance.fixed.",
      call. = FALSE
    )
  }

  dr.method <- match.arg(dr.method)

  n <- length(a)
  x_sop <- as_matrix.cir(x_sop, n = n, name = "x_sop")
  x_ps <- as_matrix.cir(x_ps, n = n, name = "x_ps")
  x_cen <- as_matrix.cir(x_cen, n = n, name = "x_cen")
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  p_rho <- fun_tx_dim.cir(fun_tx)
  p_tau <- ncol(x_sop)

  outcome.estimated <- TRUE
  if (is.null(rho.nuis) || is.null(tau.nuis)) {
    if (is.null(mle.fit)) {
      mle.fit <- MLEst.cir.analytic(
        a, fun_tx, x_sop, y, delta_c,
        rho.start = rho.start, tau.start = tau.start,
        max.step = max.step, thres = thres, timepoints = timepoints,
        vcov.type = mle.vcov.type
      )
    }
    rho.nuis <- mle.fit$point.est[seq_len(p_rho)]
    tau.nuis <- mle.fit$point.est[p_rho + seq_len(p_tau)]
  } else if (is.null(mle.fit)) {
    if (!supplied.nuisance.fixed) {
      stop(
        "rho.nuis/tau.nuis were supplied without mle.fit. Supply the fit that ",
        "generated them, or set supplied.nuisance.fixed = TRUE.",
        call. = FALSE
      )
    }
    outcome.estimated <- FALSE
  }

  ps.estimated <- TRUE
  if (is.null(gamma.nuis)) {
    if (is.null(ps.fit)) ps.fit <- fit_ps.cir(a, x_ps)
    gamma.nuis <- ps.fit$point.est
  } else if (is.null(ps.fit)) {
    if (!supplied.nuisance.fixed) {
      stop(
        "gamma.nuis was supplied without ps.fit. Supply the fit that generated ",
        "it, or set supplied.nuisance.fixed = TRUE.",
        call. = FALSE
      )
    }
    ps.estimated <- FALSE
  }

  cen.estimated <- TRUE
  if (is.null(eta.nuis)) {
    if (is.null(cen.fit)) {
      cen.fit <- fit_censoring.cir.analytic(
        y, delta_c, x_cen, eta.start = eta.start,
        max.step = max.step, thres = thres,
        time.dep = time.dep, timepoints = timepoints
      )
    }
    eta.nuis <- cen.fit$point.est
  } else if (is.null(cen.fit)) {
    if (!supplied.nuisance.fixed) {
      stop(
        "eta.nuis was supplied without cen.fit. Supply the fit that generated ",
        "it, or set supplied.nuisance.fixed = TRUE.",
        call. = FALSE
      )
    }
    cen.estimated <- FALSE
  }

  if (identical(dr.method, "checked")) {
    dr <- dr.estimate.cir.checked(
      a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
      rho.nuis, tau.nuis, gamma.nuis, eta.nuis,
      time.dep = time.dep,
      rho.start = if (is.null(rho.start)) rho.nuis else rho.start,
      max.step = max.step, thres = thres,
      timepoints = timepoints,
      cdf.tol = cdf.tol, score.tol = score.tol
    )
  } else {
    dr <- dr.estimate.cir(
      a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
      rho.nuis, tau.nuis, gamma.nuis, eta.nuis,
      time.dep = time.dep,
      rho.start = if (is.null(rho.start)) rho.nuis else rho.start,
      max.step = max.step, thres = thres,
      timepoints = timepoints
    )
  }

  var <- var.DR.sandwich.analytic.cir(
    dr.point = dr$point.esti,
    a = a, fun_tx = fun_tx, x_sop = x_sop, x_ps = x_ps, x_cen = x_cen,
    y = y, delta_c = delta_c,
    rho.nuis = rho.nuis, tau.nuis = tau.nuis,
    gamma.nuis = gamma.nuis, eta.nuis = eta.nuis,
    mle = mle.fit, ps = ps.fit, cen = cen.fit,
    outcome.estimated = outcome.estimated,
    ps.estimated = ps.estimated,
    cen.estimated = cen.estimated,
    time.dep = time.dep, timepoints = timepoints,
    center.meat = center.meat, hc1 = hc1, check_cdf = check_cdf,
    condition.warn = condition.warn
  )

  se <- var$sd.est
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  nm <- paste0("rho", seq_len(p_rho) - 1L)
  names(dr$point.esti) <- names(se) <- nm

  list(
    dr.point = dr$point.esti,
    sd.est = se,
    cov.mat = var$cov.mat,
    conf.lower = dr$point.esti - z * se,
    conf.upper = dr$point.esti + z * se,
    p.value = 2 * stats::pnorm(-abs(dr$point.esti / se)),
    convergence = dr$convergence,
    dr.method = dr.method,
    estimating_equation_max_abs = max(abs(colMeans(dr$S_mat))),
    dr_optim_convergence = cir_or_null.cir(dr$optim.convergence, dr$convergence),
    cdf.tol = cdf.tol,
    score.tol = score.tol,
    S_mat = dr$S_mat,
    mle = mle.fit,
    ps = ps.fit,
    cen = cen.fit,
    sandwich = var,
    nuisance_estimated = c(
      outcome = outcome.estimated,
      propensity = ps.estimated,
      censoring = cen.estimated
    )
  )
}
