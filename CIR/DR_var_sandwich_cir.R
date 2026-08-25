################################################################################
######################## CIR DR sandwich variance ##############################
################################################################################

cir_sum_score_for_params <- function(pars, dr.point, a, fun_tx, x_sop, x_ps, x_cen,
                                     y, delta_c, dims, time.dep, timepoints) {
  idx <- 0L
  rho.nuis <- pars[idx + seq_len(dims$p_rho)]
  idx <- idx + dims$p_rho
  tau.nuis <- pars[idx + seq_len(dims$p_tau)]
  idx <- idx + dims$p_tau
  gamma.nuis <- pars[idx + seq_len(dims$p_gamma)]
  idx <- idx + dims$p_gamma
  eta.nuis <- pars[idx + seq_len(dims$p_eta)]
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  F_nuis <- estimate_cif0.cir(rho.nuis, tau.nuis, fun_tx, x_sop, timepoints)
  e_vec <- expit.cir(drop(x_ps %*% gamma.nuis))
  G_mat <- predict_censoring_gG.cir(y, x_cen, eta.nuis, time.dep, timepoints)$G
  colSums(cir_score_matrix(dr.point, a, fun_tx, F_nuis$F0, e_vec, G_mat,
                           y, delta_c, timepoints, rho.weight = rho.nuis))
}

cir_sum_score_for_rho <- function(rho, rho.nuis, tau.nuis, gamma.nuis, eta.nuis,
                                  a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
                                  time.dep, timepoints) {
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  F_nuis <- estimate_cif0.cir(rho.nuis, tau.nuis, fun_tx, x_sop, timepoints)
  e_vec <- expit.cir(drop(x_ps %*% gamma.nuis))
  G_mat <- predict_censoring_gG.cir(y, x_cen, eta.nuis, time.dep, timepoints)$G
  colSums(cir_score_matrix(rho, a, fun_tx, F_nuis$F0, e_vec, G_mat,
                           y, delta_c, timepoints, rho.weight = rho.nuis))
}

var.DR.sandwich.cir <- function(dr.point, a, fun_tx, x_sop, x_ps, x_cen,
                                y, delta_c, rho.nuis, tau.nuis, gamma.nuis,
                                eta.nuis, S_mat_dr, mle, ps, cen,
                                time.dep = FALSE, timepoints = NULL) {
  timepoints <- resolve_timepoints.cir(y, fun_tx, timepoints)
  dims <- list(
    p_rho = length(rho.nuis),
    p_tau = length(tau.nuis),
    p_gamma = length(gamma.nuis),
    p_eta = length(eta.nuis)
  )
  nuisance <- c(rho.nuis, tau.nuis, gamma.nuis, eta.nuis)
  dU.by.nuis <- numDeriv::jacobian(function(pars) {
    cir_sum_score_for_params(pars, dr.point, a, fun_tx, x_sop, x_ps, x_cen,
                             y, delta_c, dims, time.dep, timepoints)
  }, nuisance)
  dU.by.rho <- numDeriv::jacobian(function(rho) {
    cir_sum_score_for_rho(rho, rho.nuis, tau.nuis, gamma.nuis, eta.nuis,
                          a, fun_tx, x_sop, x_ps, x_cen, y, delta_c,
                          time.dep, timepoints)
  }, dr.point)

  i1 <- seq_len(dims$p_rho + dims$p_tau)
  i2 <- max(i1) + seq_len(dims$p_gamma)
  i3 <- max(i2) + seq_len(dims$p_eta)

  U_tilde <- S_mat_dr
  U_tilde <- U_tilde + t(dU.by.nuis[, i1, drop = FALSE] %*% mle$cov %*% t(mle$S_mat))
  U_tilde <- U_tilde + t(dU.by.nuis[, i2, drop = FALSE] %*% ps$cov %*% t(ps$S_mat))
  U_tilde <- U_tilde + t(dU.by.nuis[, i3, drop = FALSE] %*% cen$cov %*% t(cen$S_mat))

  tau.inv <- -safe_solve.cir(dU.by.rho)
  cov.mat <- tau.inv %*% (t(U_tilde) %*% U_tilde) %*% t(tau.inv)
  sd.est <- sqrt(diag(cov.mat))
  list(
    sd.est = sd.est,
    cov.mat = cov.mat,
    U_tilde = U_tilde,
    dU.by.rho = dU.by.rho,
    dU.by.nuis = dU.by.nuis
  )
}
