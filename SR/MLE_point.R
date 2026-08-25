################################################################################
####################### Point Estimation Via MLE ###############################
################################################################################

max.likelihood = function(a, x_a = NULL, x_b, y, delta_c, alpha.start, beta.start,
                          max.step, thres, fun_tx = NULL){
  ## ---------------------------------------------------------------------------
  ## The name of the function: max.likelihood
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to estimate (alpha, beta) via maximum likelihood without
  ## censoring.
  ## ---------------------------------------------------------------------------
  ## Input:
  ## a: The exposure vector. Should only take values 0 or 1.
  ## x_a, x_b: covariate matrices (all the elements in the first column are 1).
  ## y: The observed time vector.
  ## delta_c: event indicator vector.
  ## alpha.start: Starting values for the parameters in the target model. ((p+1)-vector)
  ## beta.start: Starting values for the parameters in the nuisance model. ((p+1)-vector)
  ## max.step: The maximum number of iterations to be passed into the optim function.
  ## thres: Threshold for judging covergence.
  ## ---------------------------------------------------------------------------
  ## Output:
  ## opt: a list containing estimated (alpha, beta).
  ## ---------------------------------------------------------------------------


  startpars = c(alpha.start, beta.start)
  timepoints = sort(unique(c(0,y)))
  max_int_diff_time = diff(timepoints)
  min_abs_log_diff_time = min(abs(log(diff(timepoints))))
  fun_tx = .prepare.fun_tx(x_a = x_a, fun_tx = fun_tx, timepoints = timepoints)
  pa = ncol(fun_tx[[1]])

  if (is.vector(x_b)){
    pb = 1
  }else{
    pb = dim(x_b)[2]
  }

  max_abs_delta_fun = max(abs(unlist(lapply(1:(length(timepoints) - 1), function(j) {
    fun_tx[[j + 1]] - fun_tx[[j]]
  }))), na.rm = TRUE)
  if (!is.finite(max_abs_delta_fun) || max_abs_delta_fun == 0) {
    max_abs_delta_fun = 1
  }

  ### negative log likelihood function
  neg.log.likelihood = function(pars){
    n = length(a)

    alpha = pars[1:pa]
    beta = pars[(pa + 1):(pa + pb)]

    p0p1 = getProb.sequent(alpha, beta, x_b = x_b, timepoints = timepoints, fun_tx = fun_tx) # length(p0p1) = m

    NeLogL = 0
    for (i in 1:n){
      a_i = a[i]
      y_i = y[i]
      delta_c_i = delta_c[i]
      index_timepoints = which(timepoints == y_i)

      p_ai = NULL
      for(j in 1:(index_timepoints - 1)){
        p_ai = c(p_ai, p0p1[[j]][i,a_i+1])
      }

      if (length(p_ai) > 1){
        NeLogL = NeLogL - delta_c_i * log(prod(p_ai[1:(index_timepoints - 2)]) - prod(p_ai)) -
          (1 - delta_c_i) * log(prod(p_ai))
      }else{
        NeLogL = NeLogL - delta_c_i * log(1 - p_ai) - (1 - delta_c_i) * log(p_ai)
      }

    }

    return(NeLogL)
  }

  neg.log.likelihood.i = function(pars){
    n = length(a)

    alpha = pars[1:pa]
    beta = pars[(pa + 1):(pa + pb)]

    p0p1 = getProb.sequent(alpha, beta, x_b = x_b, timepoints = timepoints, fun_tx = fun_tx) # length(p0p1) = m

    NeLogL = c()
    for (i in 1:n){
      a_i = a[i]
      y_i = y[i]
      delta_c_i = delta_c[i]
      index_timepoints = which(timepoints == y_i)

      p_ai = NULL
      for(j in 1:(index_timepoints - 1)){
        p_ai = c(p_ai, p0p1[[j]][i,a_i+1])
      }

      if (length(p_ai) > 1){
        NeLogL = c(NeLogL, - delta_c_i * log(prod(p_ai[1:(index_timepoints - 2)]) - prod(p_ai)) -
                     (1 - delta_c_i) * log(prod(p_ai)))
      }else{
        NeLogL = c(NeLogL, - delta_c_i * log(1 - p_ai) - (1 - delta_c_i) * log(p_ai))
      }

    }

    return(NeLogL)
  }

  neg.log.likelihood.alpha = function(alpha){
    n = length(a)

    p0p1 = getProb.sequent(alpha, beta, x_b = x_b, timepoints = timepoints, fun_tx = fun_tx)

    NeLogL = 0
    for (i in 1:n){
      a_i = a[i]
      y_i = y[i]
      delta_c_i = delta_c[i]
      index_timepoints = which(timepoints == y_i)

      p_ai = NULL
      for(j in 1:(index_timepoints - 1)){
        p_ai = c(p_ai, p0p1[[j]][i,a_i+1])
      }

      if (length(p_ai) > 1){
        NeLogL = NeLogL - delta_c_i * log(prod(p_ai[1:(index_timepoints - 2)]) - prod(p_ai)) -
          (1 - delta_c_i) * log(prod(p_ai))
      }else{
        NeLogL = NeLogL - delta_c_i * log(1 - p_ai) - (1 - delta_c_i) * log(p_ai)
      }
      # print(c(i,NeLogL))
    }

    return(NeLogL)
  }


  neg.log.likelihood.beta = function(beta){
    n = length(a)

    p0p1 = getProb.sequent(alpha, beta, x_b = x_b, timepoints = timepoints, fun_tx = fun_tx)

    NeLogL = 0
    for (i in 1:n){
      a_i = a[i]
      y_i = y[i]
      delta_c_i = delta_c[i]
      index_timepoints = which(timepoints == y_i)

      p_ai = NULL
      for(j in 1:(index_timepoints - 1)){
        p_ai = c(p_ai, p0p1[[j]][i,a_i+1])
      }

      if (length(p_ai) > 1){
        NeLogL = NeLogL - delta_c_i * log(prod(p_ai[1:(index_timepoints - 2)]) - prod(p_ai)) -
          (1 - delta_c_i) * log(prod(p_ai))
      }else{
        NeLogL = NeLogL - delta_c_i * log(1 - p_ai) - (1 - delta_c_i) * log(p_ai)
      }

    }

    return(NeLogL)
  }


  ### Optimization
  Diff = function(z,y) sum((z - y)^2)/sum(z^2 + thres)
  alpha = alpha.start
  beta = beta.start
  diff = thres + 1
  step = 0

  while(diff > thres & step < max.step){
    step = step + 1
    if (length(alpha.start) > 1){
      opt1 = stats::optim(alpha, neg.log.likelihood.alpha, control = list(maxit = max(100, max.step/10)))
    }else{
      opt1 = stats::optim(alpha, neg.log.likelihood.alpha, control = list(maxit = max(100, max.step/10)), method = "Brent", lower = -12/max_abs_delta_fun, upper = 12/max_abs_delta_fun)
    }
    diff1 = Diff(opt1$par, alpha)
    alpha = opt1$par
    if (length(beta.start) > 1){
      opt2 = stats::optim(beta, neg.log.likelihood.beta, control = list(maxit = max(100, max.step/10)))
    }else{
      opt2 = stats::optim(beta, neg.log.likelihood.beta, control = list(maxit = max(100, max.step/10)), method = "Brent", lower = min_abs_log_diff_time - 12, upper = 12 - min_abs_log_diff_time)
    }
    diff = max(diff1, Diff(opt2$par, beta))
    beta = opt2$par
  }

  # Hessian and score at the final estimate
  point <- c(alpha, beta)
  H <- numDeriv::hessian(func = neg.log.likelihood, x = point)
  cov.optim <- solve(H)
  S_mat <- -numDeriv::jacobian(func = neg.log.likelihood.i, x = point)  # rows = obs, cols = params; n by p

  opt = list(par = point, 
             convergence = (diff <= thres && opt1$convergence == 0 && 
                              opt2$convergence == 0),
             value = neg.log.likelihood(point), step = step,
             cov.optim = cov.optim,
             S_mat = S_mat)

  

  return(opt)
}



