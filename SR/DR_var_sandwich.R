## Sandwich estimator for variance
## all sub-parts are including n, so no need to scale by n afterwards

hessian.DR <- function(a, x_a = NULL, x_b_nuis, x_c_nuis, y, delta_c,
                   alpha.dr, alpha.nuis, beta.nuis, gamma.nuis, fun_tx = NULL){
  ### U = w * R * (A - e(X))
  ### w: alpha.nuis, beta.nuis, gamma.nuis
  ### R: alpha.dr, alpha.nuis, beta.nuis
  ### e(X): gamma.nuis
  ### cov.lop: -E^{-1}(d S / d(alpha,beta))
  timepoints = sort(unique(c(0,y)))
  m <- length(timepoints) - 1
  n <- length(a)
  fun_tx <- .prepare.fun_tx(x_a = x_a, fun_tx = fun_tx, timepoints = timepoints)
  p <- ncol(fun_tx[[1]])
  pa <- ncol(fun_tx[[1]])
  pb <- dim(x_b_nuis)[2]
  pc <- dim(x_c_nuis)[2]


  # cov.psc: -E^{-1}(d S / d gamma)
  # and score_mat_psc
  e_vec <- as.vector(expit(x_c_nuis %*% gamma.nuis))
  S_mat_psc <- x_c_nuis * (a - e_vec) # n by pc
  dpscore.by.dgamma = x_c_nuis * e_vec * (1 - e_vec)  # n by pc
  hessian.psc <- -t(x_c_nuis) %*% dpscore.by.dgamma  # pc by pc
  cov.psc <- -solve(hessian.psc)

  Tau <- diag(m)
  if (m > 1) {
    for(k in 1:(m-1)) Tau[k, k+1] <- -1
  }

  # E(d U / d all)
  ## efficient score total summation
  effscore.sum <- function(pars){
    alpha <- pars[1:pa]
    alpha.nuis <- pars[(pa+1):(pa+pa)]
    beta.nuis <- pars[(2*pa+1):(2*pa+pb)]
    gamma.nuis <- pars[(2*pa+pb+1):(2*pa+pb+pc)]

    p0p1 <- getProb.sequent(alpha.nuis, beta.nuis, x_b = x_b_nuis,
                            timepoints = timepoints, fun_tx = fun_tx)

    out <- numeric(p)

    for (i in seq_len(n)) {
      e_i   <- expit(drop(x_c_nuis[i, ] %*% gamma.nuis))
      A_i   <- a[i]
      Y_i   <- y[i]
      dlt_i <- delta_c[i]

      # F_i[k,] = f_sr(t_k; X_i) for k = 1,...,m
      F_i <- do.call(rbind, lapply(2:(m+1), function(k) fun_tx[[k]][i,, drop = FALSE]))

      # streaming accumulators
      S0_prev <- 1
      S1_prev <- 1
      H_prev  <- 1
      u_vec   <- numeric(m)

      for (k in seq_len(m)) {
        t_k1 <- timepoints[k + 1]

        # update survival factors using per-interval hazards
        p0 <- p0p1[[k]][i, 1]
        p1 <- p0p1[[k]][i, 2]
        S0_k <- S0_prev * p0
        S1_k <- S1_prev * p1

        # q0_k, q1_k, dx_k (handled separately for k=1 vs k>=2)
        if (k == 1) {
          denom0 <- 1 - S0_k
          q0_k   <- S0_k / denom0
          q1_k   <- S1_k / (1 - S1_k)
          dx_k   <- S0_k / (e_i * denom0)
        } else {
          r0   <- S0_k / S0_prev
          r1   <- S1_k / S1_prev
          q0_k <- S0_k / (1 - r0)
          q1_k <- S1_k / (1 - r1)
          dx_k <- S0_k / (e_i * (1 - r0))
        }

        # mixing ratio
        num <- q1_k * e_i
        den <- num + q0_k * (1 - e_i)
        ratioA <- if (den > 0) num / den else 0

        # r_k (risk and censor pieces); keep original logic for risk set
        eta    <- A_i * sum(alpha * F_i[k, ])
        riskset <- all(Y_i != timepoints[1:(k + 1)]) # TRUE iff Y_i > t[k]
        H_k    <- exp(-eta) * as.numeric(riskset)
        C_k    <- exp(-eta) * (1 - dlt_i) * as.numeric(Y_i == t_k1)

        r_k <- H_k / S0_k + C_k / S0_k - H_prev / S0_prev

        # u_k for this interval
        u_vec[k] <- dx_k * ratioA * (A_i - e_i) * r_k

        # roll forward
        S0_prev <- S0_k
        S1_prev <- S1_k
        H_prev  <- H_k
      }

      out <- out + as.vector(t(F_i) %*% Tau %*% u_vec)
    }

    out
  }

  cov.effscore.large <- jacobian(func = effscore.sum,
                                 x     = c(alpha.dr, alpha.nuis, beta.nuis, gamma.nuis),
                                 method = "Richardson")


  ##
  dU.by.dalpha <- cov.effscore.large[,1:pa]
  dU.by.dab <- cov.effscore.large[,(pa+1):(2*pa+pb)]
  dU.by.dgamma <- cov.effscore.large[,(2*pa+pb+1):(2*pa+pb+pc)]

  return(list(S_mat_psc=S_mat_psc, cov.psc=cov.psc,
              dU.by.dalpha=dU.by.dalpha,
              dU.by.dab=dU.by.dab,
              dU.by.dgamma=dU.by.dgamma))

}



var.DR.sandwich <- function(dr.point, S_mat_psc, S_mat_lop, S_mat_dr,
                            cov.psc, cov.lop,
                            dU.by.dalpha, dU.by.dab, dU.by.dgamma){
  U_tilde <- S_mat_dr + t(dU.by.dgamma %*% cov.psc %*% t(S_mat_psc)) +
    t(dU.by.dab %*% cov.lop %*% t(S_mat_lop))  # n by p

  tau.inv <- -solve(dU.by.dalpha)  # p by p

  # sandwich variance estimation
  cov.mat <- tau.inv %*% (t(U_tilde) %*% U_tilde) %*% t(tau.inv)
  sd.est <- sqrt(diag(cov.mat))


  conf.lower <- dr.point + stats::qnorm(0.025) * sd.est
  conf.upper <- dr.point + stats::qnorm(0.975) * sd.est
  p.temp <- stats::pnorm(dr.point/sd.est, 0, 1)
  p.value <- 2 * pmin(p.temp, 1 - p.temp)

  return(list(sd.est = sd.est, cov.mat = cov.mat,
              conf.lower = conf.lower, conf.upper = conf.upper, p.value = p.value))



}
