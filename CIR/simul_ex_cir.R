################################################################################
#################### Simulation example for CIR chapter 4 ######################
################################################################################

Simul.Ex.cir <- function(rho.tr, tau.tr, gamma.tr = c(0.1, -0.5),
                         timepoints, n, sample.seed = NULL,
                         eta.tr = c(-2.5, 0.3)) {
  if (!is.null(sample.seed)) set.seed(sample.seed)
  m <- length(timepoints)
  p <- length(rho.tr) - 1L
  x_pre <- matrix(stats::runif(n * p, -2, 2), n, p)
  x <- cbind(1, x_pre)

  fun_tx <- vector("list", m + 1L)
  for (k in seq_len(m)) {
    fun_tx[[k]] <- (1 / timepoints[k]) * x
  }
  fun_tx[[m + 1L]] <- matrix(0, nrow = n, ncol = length(rho.tr))
  attr(fun_tx, "timepoints") <- timepoints

  p0p1 <- getProb.sequent.cir(rho.tr, tau.tr, fun_tx, x, timepoints)
  ps_true <- expit.cir(drop(x %*% gamma.tr))
  a <- stats::rbinom(n, 1, ps_true)

  seq_prob <- matrix(NA_real_, n, m)
  for (k in seq_len(m)) {
    seq_prob[a == 0, k] <- p0p1[[k]][a == 0, 1]
    seq_prob[a == 1, k] <- p0p1[[k]][a == 1, 2]
  }

  probs <- matrix(0, n, m + 1L)
  F_tail <- matrix(1, n, m + 1L)
  for (k in m:1) {
    F_tail[, k] <- F_tail[, k + 1L] * seq_prob[, k]
  }
  probs[, 1] <- F_tail[, 1]
  if (m > 1L) {
    for (k in 2:m) {
      probs[, k] <- F_tail[, k] - F_tail[, k - 1L]
    }
  }
  probs[, m + 1L] <- 1 - F_tail[, m]
  probs <- probs / rowSums(probs)

  t_event <- numeric(n)
  for (i in seq_len(n)) {
    idx <- sample.int(m + 1L, size = 1L, prob = probs[i, ])
    t_event[i] <- if (idx <= m) timepoints[idx] else Inf
  }

  c_time <- rep(max(timepoints), n)
  at_risk <- rep(TRUE, n)
  for (k in seq_len(m)) {
    idx <- which(at_risk)
    if (length(idx) == 0L) break
    linpred <- drop(x[idx, , drop = FALSE] %*% eta.tr)
    p_ck <- expit.cir(linpred)
    cens_now <- idx[stats::rbinom(length(idx), 1, p_ck) == 1L]
    if (length(cens_now) > 0L) {
      c_time[cens_now] <- timepoints[k]
      at_risk[cens_now] <- FALSE
    }
  }

  delta_c <- as.integer(t_event <= c_time)
  y <- ifelse(delta_c == 1L, t_event, c_time)
  y[is.infinite(y)] <- max(timepoints)
  list(a = a, y = y, delta_c = delta_c, fun_tx = fun_tx, x = x,
       ps_true = ps_true, t_event = t_event, c_time = c_time)
}
