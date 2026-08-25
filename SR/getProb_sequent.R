################################################################################
############################## getProb.sequent #################################
################################################################################

.as.design.matrix = function(x, name = "design") {
  if (is.null(x)) {
    stop(paste0(name, " must not be NULL."))
  }
  if (is.null(dim(x))) {
    matrix(x, ncol = 1)
  } else {
    as.matrix(x)
  }
}

.prepare.fun_tx = function(x_a = NULL, fun_tx = NULL, timepoints) {
  if (length(timepoints) < 2) {
    stop("timepoints must contain at least two values.")
  }

  if (!is.null(fun_tx)) {
    if (!is.list(fun_tx)) {
      stop("fun_tx must be a list with one design vector/matrix per time point.")
    }
    if (length(fun_tx) != length(timepoints)) {
      stop("length(fun_tx) must equal length(timepoints).")
    }

    fun_tx = lapply(seq_along(fun_tx), function(k) {
      .as.design.matrix(fun_tx[[k]], paste0("fun_tx[[", k, "]]"))
    })

    n = nrow(fun_tx[[1]])
    p = ncol(fun_tx[[1]])
    for (k in seq_along(fun_tx)) {
      if (nrow(fun_tx[[k]]) != n || ncol(fun_tx[[k]]) != p) {
        stop("All elements of fun_tx must have the same dimensions.")
      }
    }
    return(fun_tx)
  }

  x_a = .as.design.matrix(x_a, "x_a")
  lapply(timepoints, function(tt) x_a * tt)
}

getProb.sequent = function(alpha, beta, x_a = NULL, x_b, timepoints, fun_tx = NULL){
  ## ---------------------------------------------------------------------------
  ## The name of the function: getProb.sequent
  ## ---------------------------------------------------------------------------
  ## Description:
  ## The function to calculate sequential survival functions.
  ## ---------------------------------------------------------------------------
  ## Input:
  ## alpha: parameters in the target model.
  ## beta: parameters in the nuisance model.
  ## x_a, x_b: covariate matrices in two models (all the elements in the first column are 1).
  ## fun_tx: optional list whose kth element is f_sr(t_k; X).
  ## timepoints: recorded timepoints: (t_0, t_1, ..., t_m) where t_0 = 0.
  ## ---------------------------------------------------------------------------
  ## Output:
  ## p0p1: a list storing the sequential survival functions.
  ## ---------------------------------------------------------------------------

  ### create the list to store the sequential survival functions
  # n = nrow(x_a)
  m = length(timepoints) - 1

  p0p1 = list()
  fun_tx = .prepare.fun_tx(x_a = x_a, fun_tx = fun_tx, timepoints = timepoints)

  if (is.matrix(x_b)){
    x_b_T_beta = x_b%*%beta
  }else{
    x_b_T_beta = x_b*beta
  }

  for (j in 0:(m - 1)){
    delta_fun_tx = fun_tx[[j+2]] - fun_tx[[j+1]]
    p0p1[[j+1]] = brm::getProbRR(logrr = delta_fun_tx%*%alpha,
                                 logop = (x_b_T_beta) - log(timepoints[j+2] - timepoints[j+1]))
  }
  return(p0p1)
}

