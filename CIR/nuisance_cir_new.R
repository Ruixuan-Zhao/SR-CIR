################################################################################
###################### Propensity and censoring nuisance #######################
################################################################################

fit_ps.cir <- function(a, x_ps, conf.level = 0.95) {
  x_ps <- as_matrix.cir(x_ps, n = length(a), name = "x_ps")
  fit <- stats::glm.fit(x = x_ps, y = a, family = stats::binomial(link = "logit"))
  gamma <- as.numeric(fit$coefficients)
  e <- expit.cir(drop(x_ps %*% gamma))
  S_mat <- x_ps * as.numeric(a - e)
  hessian <- -t(x_ps) %*% (x_ps * as.numeric(e * (1 - e)))
  cov <- -safe_solve.cir(hessian)
  se <- sqrt(diag(cov))
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  list(
    point.est = gamma,
    se.est = se,
    cov = cov,
    S_mat = S_mat,
    conf.lower = gamma - z * se,
    conf.upper = gamma + z * se,
    fitted = e,
    convergence = fit$converged
  )
}

censoring_loglik_i.cir <- function(eta, y, delta_c, x_cen, timepoints,
                                   time.dep = FALSE) {
  n <- length(y)
  m <- length(timepoints)
  x_cen <- as_matrix.cir(x_cen, n = n, name = "x_cen")
  lambda <- censoring_hazard.cir(x_cen, eta, m, time.dep)
  ll <- numeric(n)
  for (i in seq_len(n)) {
    k <- match(y[i], timepoints)
    if (m == 1L) {
      ll[i] <- 0
    } else if (k == m) {
      ll[i] <- sum(log(1 - lambda[i, seq_len(m - 1L)]))
    } else {
      before <- if (k > 1L) sum(log(1 - lambda[i, seq_len(k - 1L)])) else 0
      if (delta_c[i] == 0) {
        ll[i] <- before + log(lambda[i, k])
      } else {
        ll[i] <- before
      }
    }
  }
  ll
}

censoring_hazard.cir <- function(x_cen, eta, m, time.dep = FALSE) {
  n <- nrow(x_cen)
  if (m <= 1L) {
    return(matrix(numeric(0), nrow = n, ncol = 0))
  }
  if (!time.dep) {
    lambda_i <- expit.cir(drop(x_cen %*% eta))
    matrix(rep(lambda_i, times = m - 1L), nrow = n, ncol = m - 1L)
  } else {
    p_cov <- ncol(x_cen) - 1L
    alpha <- eta[seq_len(m - 1L)]
    beta <- if (p_cov > 0L) eta[(m):(m - 1L + p_cov)] else numeric(0)
    xb <- if (p_cov > 0L) drop(x_cen[, -1, drop = FALSE] %*% beta) else rep(0, n)
    expit.cir(matrix(xb, nrow = n, ncol = m - 1L) +
                              matrix(rep(alpha, each = n), nrow = n))
  }
}

fit_censoring.cir <- function(y, delta_c, x_cen, eta.start = NULL,
                              max.step = 1000, thres = 1e-6,
                              time.dep = FALSE, timepoints = NULL) {
  x_cen <- as_matrix.cir(x_cen, n = length(y), name = "x_cen")
  timepoints <- resolve_timepoints.cir(y, timepoints = timepoints)
  m <- length(timepoints)
  if (is.null(eta.start)) {
    eta.start <- if (time.dep) rep(0, (m - 1L) + ncol(x_cen) - 1L) else rep(0, ncol(x_cen))
  }
  neg_loglik <- function(eta) {
    -sum(censoring_loglik_i.cir(eta, y, delta_c, x_cen, timepoints, time.dep))
  }
  opt <- stats::optim(
    eta.start, neg_loglik,
    control = list(maxit = max.step, reltol = thres),
    hessian = TRUE
  )
  eta <- opt$par
  S_mat <- numDeriv::jacobian(function(pars) {
    censoring_loglik_i.cir(pars, y, delta_c, x_cen, timepoints, time.dep)
  }, eta)
  cov <- safe_solve.cir(opt$hessian)
  list(
    point.est = eta,
    cov = cov,
    S_mat = S_mat,
    convergence = opt$convergence,
    value = opt$value,
    timepoints = timepoints,
    time.dep = time.dep
  )
}

predict_censoring_gG.cir <- function(y, x_cen, eta, time.dep = FALSE,
                                     timepoints = NULL) {
  x_cen <- as_matrix.cir(x_cen, n = length(y), name = "x_cen")
  timepoints <- resolve_timepoints.cir(y, timepoints = timepoints)
  m <- length(timepoints)
  n <- nrow(x_cen)
  if (m == 1L) {
    return(list(g = matrix(1, n, 1), G = matrix(0, n, 1)))
  }
  lambda <- censoring_hazard.cir(x_cen, eta, m, time.dep)
  S_after <- t(apply(1 - lambda, 1, cumprod))
  S_before <- if (m - 1L == 1L) {
    matrix(1, nrow = n, ncol = 1)
  } else {
    cbind(1, S_after[, seq_len(m - 2L), drop = FALSE])
  }
  g_early <- lambda * S_before
  g_last <- S_after[, m - 1L]
  g <- cbind(g_early, g_last)
  G <- matrix(0, nrow = n, ncol = m)
  if (m > 1L) {
    cg <- t(apply(g, 1, cumsum))
    G[, 2:m] <- cg[, seq_len(m - 1L), drop = FALSE]
  }
  list(g = g, G = G)
}
