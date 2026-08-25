################################################################################
######################## Backward sequential CDF ratios ########################
################################################################################

getProb.sequent.cir <- function(rho, tau, fun_tx, x_sop, timepoints) {
  
  timepoints <- sort(unique(timepoints))
  m <- length(timepoints)
  x_sop <- as_matrix.cir(x_sop, name = "x_sop")
  n <- nrow(x_sop)
  
  f_rho <- vector("list", m + 1)
  for (j in seq_len(m + 1)) {
    f_rho[[j]] <- drop(fun_tx_matrix.cir(fun_tx, j, n = n) %*% rho)
  }
  x_tau <- drop(x_sop %*% tau)

  out <- vector("list", m)
  for (j in seq_len(m - 1)) {
    out[[j]] <- brm::getProbRR(
      logrr = f_rho[[j]] - f_rho[[j + 1]],
      logop = x_tau - log(1 / timepoints[j] - 1 / timepoints[j + 1])
    )
    out[[j]] <- as.matrix(out[[j]])
  }
  out[[m]] <- brm::getProbRR(
    logrr = f_rho[[m]] - f_rho[[m + 1]],
    logop = x_tau - log(1 / timepoints[m])
  )
  out[[m]] <- as.matrix(out[[m]])
  out
}

cif_from_sequential.cir <- function(p0p1) {
  m <- length(p0p1)
  F0F1 <- vector("list", m)
  F0F1[[m]] <- as.matrix(p0p1[[m]])
  if (m > 1) {
    for (j in (m - 1):1) {
      F0F1[[j]] <- F0F1[[j + 1]] * as.matrix(p0p1[[j]])
    }
  }
  F0F1
}

estimate_cif0.cir <- function(rho.nuis, tau.nuis, fun_tx, x_sop, timepoints) {
  p0p1 <- getProb.sequent.cir(rho.nuis, tau.nuis, fun_tx, x_sop, timepoints)
  F0F1 <- cif_from_sequential.cir(p0p1)
  F0 <- do.call(cbind, lapply(F0F1, function(z) z[, 1]))
  F1 <- do.call(cbind, lapply(F0F1, function(z) z[, 2]))
  list(F0 = F0, F1 = F1, p0p1 = p0p1)
}

candidate_cif.cir <- function(rho, fun_tx, F0) {
  n <- nrow(F0)
  m <- ncol(F0)
  eta <- matrix(NA_real_, n, m)
  for (k in seq_len(m)) {
    eta[, k] <- drop(fun_tx_matrix.cir(fun_tx, k, n = n) %*% rho)
  }
  F1 <- F0 * exp(eta)
  list(F0 = F0, F1 = F1, eta = eta)
}
