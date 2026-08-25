################################################################################
################## Point Estimation Via DR estimation equation #################
################################################################################
library(brm)

expit = function(logodds) {
  1/(1 + exp(-logodds))
}

dr.estimate = function(a, x_a = NULL, x_b_nuis, x_c_nuis, y, delta_c, alpha.nuis, beta.nuis, gamma.nuis,
                          max.step, thres, alpha.start, fun_tx = NULL){

  timepoints = sort(unique(c(0,y)))
  m <- length(timepoints) - 1
  n <- length(a)
  fun_tx <- .prepare.fun_tx(x_a = x_a, fun_tx = fun_tx, timepoints = timepoints)
  p = ncol(fun_tx[[1]])


  # --- 0) compute nuisance functions
  ## p0p1
  p0p1 <- getProb.sequent(alpha.nuis, beta.nuis, x_b = x_b_nuis,
                          timepoints = timepoints, fun_tx = fun_tx)
  ## S_0 and S_1
  S0S1 <- list()
  S0S1[[1]] <- p0p1[[1]]
  if (m > 1) {
    for(j in 2:m){
      S0S1[[j]] <- S0S1[[j-1]] * p0p1[[j]]
    }
  }

  ## Propensity score e(X)
  e_vec <- as.vector(expit(x_c_nuis %*% gamma.nuis))


  # --- 1) build Tau: m by m with 1 on diag, -1 on super-diag -------------
  Tau <- diag(m)
  if (m > 1) {
    for(k in 1:(m-1)) Tau[k, k+1] <- -1
  }

  # --- 2) precompute D_X diagonal entries for each i,k -------------
  dx <- matrix(NA, n, m)
  dx[,1] <- S0S1[[1]][,1] / ( e_vec * (1 - S0S1[[1]][,1]) )
  if (m > 1) {
    for(k in 2:m){
      dx[,k] <- S0S1[[k]][,1] / ( e_vec * (1 - S0S1[[k]][,1]/S0S1[[k-1]][,1]) )
    }
  }

  # --- 3) precompute q1 and q0 for each k  ---
  q1_k  <- matrix(NA, n, m)
  q0_k  <- matrix(NA, n, m)
  q0_k[,1] <- S0S1[[1]][,1] / (1 - S0S1[[1]][,1])
  q1_k[,1] <- S0S1[[1]][,2] / (1 - S0S1[[1]][,2])
  if (m > 1) {
    for(k in 2:m){
      q0_k[,k] <- S0S1[[k]][,1] / (1 - S0S1[[k]][,1] / S0S1[[k-1]][,1])
      q1_k[,k] <- S0S1[[k]][,2] / (1 - S0S1[[k]][,2] / S0S1[[k-1]][,2])
    }
  }

  Eaq_x <- q1_k * e_vec
  Eq_x <- q1_k * e_vec + q0_k * (1 - e_vec)


  # --- 4) define the score-sum-of-squares objective ----------------
  effscore.objective <- function(alpha){
    # will fill a p by n matrix of individual scores
    S_mat <- matrix(0, p, n)

    for(i in 1:n){
      # 4a) assemble the m by p F_i matrix

      F_i <- do.call(rbind, lapply(2:(m+1), function(k) fun_tx[[k]][i,, drop = FALSE]))

      e_i   <- e_vec[i]
      A_i   <- a[i]
      Y_i   <- y[i]
      dlt_i <- delta_c[i]

      # 4b) build the R-vector of length m
      r_i <- numeric(m)
      H_prev <- 1
      S0_prev <- 1
      for(k in 1:m){
        # linear predictor at t[k]:
        eta_ik <- sum(alpha * F_i[k, ]) * A_i

        # H_k(alpha)
        riskset <- all(Y_i != timepoints[1:(k+1)])  # TRUE iff Y_i > t[k]
        H_k <- exp(-eta_ik) * as.numeric(riskset)

        # C_k(alpha)
        C_k <- exp(-eta_ik) * (1 - dlt_i) * as.numeric(Y_i == timepoints[k+1])

        r_i[k] <- H_k/S0S1[[k]][i,1] + C_k/S0S1[[k]][i,1] - H_prev/S0_prev

        H_prev  <- H_k
        S0_prev <- S0S1[[k]][i,1]
      }

      # 4c) form the efficient score vector v
      # ratio <- (q1_k[i,]*e_i) / ( q1_k[i,]*e_i + q0[i,]*(1-e_i) )
      ratio <- Eaq_x[i,] / Eq_x[i,]
      u     <- dx[i,] * ratio * (A_i - e_i) * r_i

      # 4d) apply F_i and Tau to get v
      v <- as.vector(t(F_i) %*% Tau %*% u)   # length p

      S_mat[,i] <- v
    }

    return(S_mat)

  }

  effscore.alphadr <- function(alpha){
    # will fill a p by n matrix of individual scores
    S_mat <- matrix(0, p, p)

    for(i in 1:n){
      # 4a) assemble the m by p F_i matrix

      F_i <- do.call(rbind, lapply(2:(m+1), function(k) fun_tx[[k]][i,, drop = FALSE]))

      e_i   <- e_vec[i]
      A_i   <- a[i]
      Y_i   <- y[i]
      dlt_i <- delta_c[i]

      # 4b) build the R-vector of length m
      r_i <- matrix(0, m, p)
      H_prev <- 1
      S0_prev <- 1
      F_i_prev <- as.numeric(fun_tx[[1]][i,, drop = TRUE])
      for(k in 1:m){
        # linear predictor at t[k]:
        eta_ik <- sum(alpha * F_i[k, ]) * A_i

        # H_k(alpha)
        riskset <- all(Y_i != timepoints[1:(k+1)])  # TRUE iff Y_i > t[k]
        H_k <- exp(-eta_ik) * as.numeric(riskset)

        # C_k(alpha)
        C_k <- exp(-eta_ik) * (1 - dlt_i) * as.numeric(Y_i == timepoints[k+1])

        # r_i[k] <- H_k/S0S1[[k]][i,1] + C_k/S0S1[[k]][i,1] - H_prev/S0_prev
        # change to derivative
        r_i[k, ] <- -H_k * A_i * F_i[k, ]/S0S1[[k]][i,1] - C_k * A_i * F_i[k, ]/S0S1[[k]][i,1] + H_prev * A_i * F_i_prev/S0_prev # p-dim vector

        H_prev  <- H_k
        S0_prev <- S0S1[[k]][i,1]
        F_i_prev <- F_i[k, ]
      }

      # 4c) form the efficient score vector v
      # ratio <- (q1_k[i,]*e_i) / ( q1_k[i,]*e_i + q0[i,]*(1-e_i) )
      ratio <- Eaq_x[i,] / Eq_x[i,]
      u     <- dx[i,] * ratio * (A_i - e_i) * r_i

      # 4d) apply F_i and Tau to get v
      #v <- as.vector(t(F_i) %*% Tau %*% u)   # length p
      v <- t(F_i) %*% Tau %*% u

      S_mat <- S_mat + v # p by p
    }

    return(S_mat)

  }

  dr.objective <- function(alpha){
    S_mat <- effscore.objective(alpha)
    # 4e) our estimating eqn is mean_i S_mat[,i] = 0;
    #     we turn that into a scalar by sum of squares:
    mscore <- rowMeans(S_mat)
    sum(mscore^2)
  }

  # --- 5) call optim to solve for alpha --------------------------
  opt <- stats::optim(
    alpha.start,
    dr.objective,
    control = list(maxit = max.step, reltol = thres)
  )


  point.esti <- opt$par
  convergence <- opt$convergence

  # --- 6) estimate efficient variance
  # score function values
  S_mat <- t(effscore.objective(point.esti)) # n by p

  return(list(point.esti = point.esti, convergence = convergence,
              S_mat = S_mat))


}
