################################################################################
############################ Simulation Example ################################
################################################################################

Simul.Ex.censoring = function(alpha.tr, beta.tr, gamma.tr, timepoints, n, 
                              sample.seed = NULL){
  
  if (!is.null(sample.seed)) set.seed(sample.seed)
  m = length(timepoints) - 1
  p = length(alpha.tr) - 1
  
  x.pre = matrix(runif(n*p, -2, 2),n,p)
  x = cbind(rep(1,n),x.pre)
  
  p0p1 = getProb.sequent(alpha.tr, beta.tr, x, x, timepoints)
  
  pscore.true = plogis(x %*% gamma.tr)
  a = rbinom(n, 1, pscore.true)
  
  pa = matrix(NA, n, m)
  for (j in 0:(m-1)){
    pa[a == 0,j+1] = p0p1[[j+1]][a == 0,1]
    pa[a == 1,j+1] = p0p1[[j+1]][a == 1,2]
  }
  
  Q.prob = matrix(0,n,m+1)
  Q.prob[,1] = 1 - pa[,1]
  for (k in 2:m){
    if(k > 2){
      Q.prob[,k] = apply(pa[,1:(k-1)],1,prod) - apply(pa[,1:k],1,prod)
    }else{
      Q.prob[,k] = pa[,1:(k-1)] - apply(pa[,1:k],1,prod)
    }
  }
  Q.prob[,m+1] = apply(pa[,1:m],1,prod)
  
  t = NULL
  for(i in 1:n){
    t = c(t, extraDistr::rcat(1, Q.prob[i,]))
  }
  
  t <- c(timepoints[-1], max(timepoints)+1)[t]
  
  
  c_raw = rexp(n, 2/timepoints[m+1])
  c = sapply(c_raw, function(x) {
    min(base::c(timepoints[timepoints >= x], Inf))
  })
  
  c[c==Inf] = max(timepoints)
  
  delta_c = rep(0,n)
  delta_c[t <= c] = 1
  
  y = c
  y[t <= c] = t[t <= c]
  
  return(list(a = a, y = y, delta_c = delta_c, x = x))
}


