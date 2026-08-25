################################################################################
################################################################################
########################## Monte-Carlo simulation ##############################
################################################################################
################################################################################
rm(list = ls())
source("getProb_sequent.R")
source("MLE_point.R")
source("CallMLE.R")
source("WrapResults.R")
source("simul_ex.R")
source("DR_point_var.R")
source("CallDR.R")
source("DR_var_sandwich.R")
library(parallel)
library(numDeriv)


# ---- parse command-line args (no extra packages needed) ----
args <- commandArgs(trailingOnly = TRUE)
kv <- strsplit(sub("^--", "", args), "=", fixed = TRUE)
arg_list <- setNames(lapply(kv, `[`, 2L), vapply(kv, `[`, "", 1L))

get_arg <- function(name, default = NULL) {
  val <- arg_list[[name]]
  if (is.null(val) || is.na(val) || identical(val, "")) default else val
}

# ---- required/optional inputs ----
n         <- as.integer(get_arg("n", 1000))
example   <- as.integer(get_arg("example", 2))
date_tag  <- format(Sys.Date(), "%m%d")
# seed_fix  <- 3047
seed_fix  <- 20250820
m         <- 500

# If filename not passed, make a sensible default
filename  <- get_arg("filename",
                     sprintf("example%d_n%d_seed%d_%s_debug", example, n, seed_fix, date_tag))

# ---- your original variables (kept for compatibility) ----
dir.create("Results", showWarnings = FALSE, recursive = TRUE)
path.save <- paste0("Results/", filename)
duration  <- Sys.time()

set.seed(seed_fix)

# ---- choose timepoints by example ----
if (example == 1L) {
  # example 1
  timepoints <- 0:10
} else if (example == 2L) {
  # example 2
  timepoints <- unique(c(0, sort(runif(10, 0, 10))))
} else {
  stop("`--example` must be 1 or 2.")
}

message(sprintf(">> Running example=%d, n=%d, m=%d, seed=%d", example, n, m, seed_fix))
message(sprintf(">> filename=%s", filename))


# ------- start from here
alpha.tr <- c(0, 1)
beta.tr  <- c(-0.5, 1)
gamma.tr <- c(0.1, -0.5)


p_alpha <- length(alpha.tr) - 1
p_beta  <- length(beta.tr)  - 1
p_gamma <- length(gamma.tr) - 1

p_total <- length(alpha.tr) + length(beta.tr) + length(gamma.tr)
settings <- c("bth","psc","opc","bad")
S <- length(settings)

## pre‑allocate just like before
point.est.ml <- array(NA, c(m, p_total, S), dimnames = list(NULL,NULL,settings))
sd.est.ml    <- array(NA, c(m, p_total, S))
lcb.est.ml   <- array(NA, c(m, p_total, S))
rcb.est.ml   <- array(NA, c(m, p_total, S))

point.est.dr <- array(NA, c(m, length(alpha.tr), S), dimnames = list(NULL,NULL,settings))

sd.est.sand.dr <- array(NA, c(m, length(alpha.tr), S))
lcb.est.sand.dr <- array(NA, c(m, length(alpha.tr), S))
rcb.est.sand.dr <- array(NA, c(m, length(alpha.tr), S))

## run all m replicates in parallel via forking
res_list <- mclapply(
  seq_len(m),
  function(j) {
    ## simulate with seed = j
    dat <- Simul.Ex.censoring(
      alpha.tr = alpha.tr, beta.tr = beta.tr, gamma.tr = gamma.tr,
      timepoints = timepoints, n = n, sample.seed = seed_fix + j
    )
    
    a       <- dat$a
    x       <- dat$x
    y       <- dat$y
    delta_c <- dat$delta_c
    
    x_b <- cbind(1, matrix(runif(n * p_beta, -2, 2), n, p_beta))
    x_c <- cbind(1, matrix(runif(n * p_gamma, -2, 2), n, p_gamma))
    
    
    ## per‐run holders
    pe_ml <- matrix(NA, p_total, S)
    sd_ml <- matrix(NA, p_total, S)
    lcb_ml<- matrix(NA, p_total, S)
    rcb_ml<- matrix(NA, p_total, S)
    
    pe_dr <- matrix(NA, length(alpha.tr), S)
    
    sd_sand_dr <- matrix(NA, length(alpha.tr), S)
    lcb_sand_dr<- matrix(NA, length(alpha.tr), S)
    rcb_sand_dr<- matrix(NA, length(alpha.tr), S)
    failures <- character()
    
    for (k in seq_along(settings)) {
      setn <- settings[k]
      if      (setn=="bth") { xb_n <- x;   xc_n <- x   }
      else if (setn=="psc") { xb_n <- x_b; xc_n <- x   }
      else if (setn=="opc") { xb_n <- x;   xc_n <- x_c }
      else                  { xb_n <- x_b; xc_n <- x_c }
      
      res <- tryCatch({
        ## 1) MLE + PS
        rml <- MLEst(
          a=a, x_a=x, x_b=xb_n,
          y=y, delta_c=delta_c,
          alpha.start=rep(0,p_alpha+1),
          beta.start =rep(0,p_beta+1),
          max.step=1000, thres=1e-6, method=1
        )
        rps <- fit_ps(a=a, x_c_nuis=xc_n)

        ## 2) Doubly‐robust
        an <- rml$point.est[1:(p_alpha+1)]
        bn <- rml$point.est[(p_alpha+2):(p_alpha+p_beta+2)]
        gn <- rps$point.est

        rdr <- DREst(
          a=a, x_a=x,
          x_b_nuis=xb_n, x_c_nuis=xc_n,
          y=y, delta_c=delta_c,
          alpha.nuis=an, beta.nuis=bn, gamma.nuis=gn,
          max.step=1000, thres=1e-6,
          alpha.start=an
        )

        ## sandwich variance estimation
        ret_hessian <- hessian.DR(a=a, x_a=x, x_b_nuis=xb_n, x_c_nuis=xc_n, y=y, delta_c=delta_c,
                               alpha.dr=rdr$dr.point, alpha.nuis=an, beta.nuis=bn, gamma.nuis=gn)
        ret_sandwich_dr <- var.DR.sandwich(dr.point=rdr$dr.point, S_mat_psc=ret_hessian$S_mat_psc,
                                           S_mat_lop=rml$S_mat, S_mat_dr=rdr$S_mat,
                                           cov.psc=ret_hessian$cov.psc, cov.lop=rml$cov,
                                           dU.by.dalpha=ret_hessian$dU.by.dalpha,
                                           dU.by.dab=ret_hessian$dU.by.dab,
                                           dU.by.dgamma=ret_hessian$dU.by.dgamma)

        list(rml = rml, rps = rps, rdr = rdr,
             sandwich = ret_sandwich_dr)
      }, error = function(e) {
        failures <<- c(failures, paste0(setn, ": ", conditionMessage(e)))
        NULL
      })

      if (!is.null(res)) {
        pe_ml[,k] <- c(res$rml$point.est,   res$rps$point.est)
        sd_ml[,k] <- c(res$rml$se.est,      res$rps$se.est)
        lcb_ml[,k]<- c(res$rml$conf.lower,  res$rps$conf.lower)
        rcb_ml[,k]<- c(res$rml$conf.upper,  res$rps$conf.upper)

        pe_dr[,k] <- res$rdr$dr.point
        sd_sand_dr[,k] <- res$sandwich$sd.est
        lcb_sand_dr[,k]<- res$sandwich$conf.lower
        rcb_sand_dr[,k]<- res$sandwich$conf.upper
      }
    }
    
    list(
      pe_ml=pe_ml, sd_ml=sd_ml, lcb_ml=lcb_ml, rcb_ml=rcb_ml,
      pe_dr=pe_dr,
      sd_sand_dr=sd_sand_dr, lcb_sand_dr=lcb_sand_dr, rcb_sand_dr=rcb_sand_dr,
      failures=failures
    )
  },
  mc.cores = max(1L, parallel::detectCores() - 1L, na.rm = TRUE)
)

for (j in seq_len(m)) {
  point.est.ml[j,,] <- res_list[[j]]$pe_ml
  sd.est.ml[j,,]    <- res_list[[j]]$sd_ml
  lcb.est.ml[j,,]   <- res_list[[j]]$lcb_ml
  rcb.est.ml[j,,]   <- res_list[[j]]$rcb_ml
  
  point.est.dr[j,,] <- res_list[[j]]$pe_dr
  
  sd.est.sand.dr[j,,]    <- res_list[[j]]$sd_sand_dr
  lcb.est.sand.dr[j,,]   <- res_list[[j]]$lcb_sand_dr
  rcb.est.sand.dr[j,,]   <- res_list[[j]]$rcb_sand_dr
  if (length(res_list[[j]]$failures) > 0L) {
    message("simulation ", j, " failures: ", paste(res_list[[j]]$failures, collapse = " | "))
  }
}


######################### summarize results #####################################

## true values for MLE params
true.ml <- c(alpha.tr, beta.tr, gamma.tr)
param_names_ml <- c(paste0("alpha",0:p_alpha),
                    paste0("beta",0:p_beta),
                    paste0("gamma",0:p_gamma))
param_names_dr <- paste0("alpha",0:p_alpha)

p_ml <- length(true.ml)
p_dr <- length(alpha.tr)

#–– allocate summary matrices ––#
bias_mean_ml  <- matrix(NA, S, p_ml, dimnames=list(settings, param_names_ml))
bias_se_ml    <- matrix(NA, S, p_ml, dimnames=list(settings, param_names_ml))
ratio_se_ml   <- matrix(NA, S, p_ml, dimnames=list(settings, param_names_ml))
coverage_ml   <- matrix(NA, S, p_ml, dimnames=list(settings, param_names_ml))

bias_mean_dr  <- matrix(NA, S, p_dr, dimnames=list(settings, param_names_dr))
bias_se_dr    <- matrix(NA, S, p_dr, dimnames=list(settings, param_names_dr))

## DR sandwich summaries
ratio_se_dr_sand <- matrix(NA, S, p_dr, dimnames=list(settings, param_names_dr))
coverage_dr_sand <- matrix(NA, S, p_dr, dimnames=list(settings, param_names_dr))

##– loop over settings to fill them ––##
for (k in seq_along(settings)) {
  ## MLE slices
  pe_ml  <- point.est.ml[,,k]   # m × p_ml
  se_ml  <- sd.est.ml[,,k]
  lcb_ml <- lcb.est.ml[,,k]
  ucb_ml <- rcb.est.ml[,,k]
  
  keep_ml <- complete.cases(pe_ml, se_ml, lcb_ml, ucb_ml)
  pe_ml  <- pe_ml[ keep_ml, , drop = FALSE ]
  se_ml  <- se_ml[ keep_ml, , drop = FALSE ]
  lcb_ml <- lcb_ml[keep_ml, , drop = FALSE ]
  ucb_ml <- ucb_ml[keep_ml, , drop = FALSE ]
  
  ## DR point estimates and sandwich inference
  pe_dr  <- point.est.dr[,,k]   # m × p_dr
  se_dr_sand  <- sd.est.sand.dr[,,k]
  lcb_dr_sand <- lcb.est.sand.dr[,,k]
  ucb_dr_sand <- rcb.est.sand.dr[,,k]
  
  keep_dr_sand <- complete.cases(pe_dr, se_dr_sand, lcb_dr_sand, ucb_dr_sand)
  pe_dr_sand_ok  <- pe_dr[ keep_dr_sand, , drop = FALSE ]
  se_dr_sand_ok  <- se_dr_sand[ keep_dr_sand, , drop = FALSE ]
  lcb_dr_sand_ok <- lcb_dr_sand[keep_dr_sand, , drop = FALSE ]
  ucb_dr_sand_ok <- ucb_dr_sand[keep_dr_sand, , drop = FALSE ]
  
  ## 1) MLE summaries
  # bias
  bias_mat_ml <- sweep(pe_ml, 2, true.ml, "-")
  bias_mean_ml[k, ] <- colMeans(bias_mat_ml)
  bias_se_ml[k, ]   <- apply(bias_mat_ml, 2, sd) / sqrt(nrow(pe_ml))
  
  # SE comparison (reported / Monte Carlo)
  mc_se_ml  <- apply(pe_ml, 2, sd)
  rep_se_ml <- colMeans(se_ml)
  ratio_se_ml[k, ] <- ifelse(mc_se_ml == 0, NA_real_, rep_se_ml / mc_se_ml)
  
  # coverage
  coverage_ml[k, ] <- sapply(seq_len(p_ml), function(j) {
    mean(lcb_ml[,j] < true.ml[j] & ucb_ml[,j] > true.ml[j])
  })
  
  ## 2) DR summaries using sandwich inference
  bias_mat_dr <- sweep(pe_dr_sand_ok, 2, alpha.tr, "-")
  bias_mean_dr[k, ] <- colMeans(bias_mat_dr)
  bias_se_dr[k, ]   <- apply(bias_mat_dr, 2, sd) / sqrt(nrow(pe_dr_sand_ok))
  
  mc_se_dr_sand  <- apply(pe_dr_sand_ok, 2, sd)
  rep_se_dr_sand <- colMeans(se_dr_sand_ok)
  ratio_se_dr_sand[k, ] <- ifelse(mc_se_dr_sand == 0, NA_real_, rep_se_dr_sand / mc_se_dr_sand)
  
  coverage_dr_sand[k, ] <- sapply(seq_len(p_dr), function(j) {
    mean(lcb_dr_sand_ok[,j] < alpha.tr[j] & ucb_dr_sand_ok[,j] > alpha.tr[j])
  })
}

## (optional) bundle all summaries for easy return/printing
summaries <- list(
  ML = list(
    bias_mean = bias_mean_ml,
    bias_se   = bias_se_ml,
    ratio_se  = ratio_se_ml,
    coverage  = coverage_ml
  ),
  DR = list(
    bias_mean = bias_mean_dr,
    bias_se   = bias_se_dr,
    ratio_se  = ratio_se_dr_sand,
    coverage  = coverage_dr_sand
  )
)


duration <- Sys.time() - duration

# save results
if (file.exists(path.save)) {
  file.remove(list.files(path.save, full.names = TRUE))
}
# Create a new directory
dir.create(path.save, recursive = TRUE)
path.image <- paste0(path.save,'/Simu.RData')
path.out <- paste0(path.save,'/simu.out')
save.image(file=path.image)

# save parameter setting and results
options(width = 100)
sink(file=path.out, append=TRUE)
date()
duration
cat('nsims', '\n')
m
cat('sample size', '\n')
n
cat("example\n")
example
cat('metrics', '\n')
summaries
cat('\n\n')
sink()

