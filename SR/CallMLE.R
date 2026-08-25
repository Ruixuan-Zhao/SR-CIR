################################################################################
#################################### MLE #######################################
################################################################################

MLEst = function(a, x_a = NULL, x_b, y, delta_c, alpha.start = NULL, beta.start = NULL,
                 max.step, thres, method, fun_tx = NULL){

  timepoints = sort(unique(c(0,y)))
  fun_tx.resolved = .prepare.fun_tx(x_a = x_a, fun_tx = fun_tx, timepoints = timepoints)
  pa = ncol(fun_tx.resolved[[1]]) - 1

  if (is.vector(x_b)){
    pb = 0
  }else{
    pb = dim(x_b)[2] - 1
  }




  ### starting values for parameter optimization
  if (is.null(alpha.start)){
    alpha.start = rep(0, pa + 1)
  }
  if (is.null(beta.start)){
    beta.start = rep(0, pb + 1)
  }

  ### point estimate
  mle = max.likelihood(a, x_a = x_a, x_b = x_b, y = y, delta_c = delta_c,
                       alpha.start = alpha.start, beta.start = beta.start,
                       max.step = max.step, thres = thres, fun_tx = fun_tx.resolved)
  point.est = mle$par
  coverged = mle$convergence

  alpha.ml = point.est[1:(pa+1)]
  beta.ml = point.est[(pa+2):(pa + pb + 2)]

  cov.matrix <- mle$cov.optim
  S_mat = mle$S_mat


  name = paste(c(rep("alpha", pa+1), rep("beta", pb+1)), c(1:(pa+1), 1:(pb+1)))
  sol = WrapResults(point.est  = point.est, cov = cov.matrix, S_mat = S_mat,
                    name = name, x_a = x_a, x_b = x_b, converged = coverged,
                    fun_tx = if (is.null(fun_tx)) NULL else fun_tx.resolved)

  return(sol)
}

