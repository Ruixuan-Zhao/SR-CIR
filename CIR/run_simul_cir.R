################################################################################
########################### CIR chapter 4 simulation ###########################
################################################################################

.libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))
boot_parse_args.cir <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L) {
    return(list())
  }
  kv <- strsplit(sub("^--", "", args), "=", fixed = TRUE)
  setNames(lapply(kv, `[`, 2L), vapply(kv, `[`, "", 1L))
}

boot_get_arg.cir <- function(arg_list, name, default = NULL) {
  val <- arg_list[[name]]
  if (is.null(val) || is.na(val) || identical(val, "")) default else val
}

infer_code_dir.cir <- function() {
  cmd_args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

boot_args <- boot_parse_args.cir()
code_dir <- boot_get_arg.cir(boot_args, "code_dir", Sys.getenv("CIR_CODE_DIR", unset = NA_character_))
if (is.na(code_dir) || !nzchar(code_dir)) {
  code_dir <- infer_code_dir.cir()
}
code_dir <- normalizePath(code_dir, winslash = "/", mustWork = TRUE)
setwd(code_dir)


source_files <- c(
  "helpers_cir.R",
  "getProb_sequent_cir.R",
  "MLE_point_cir.R",
  "CallMLE_cir.R",
  "nuisance_cir_new.R",
  "DR_var_sandwich_cir.R",
  "DR_point_cir.R",
  "simul_ex_cir.R",
  "CIR_analytic_patch.R"
)
invisible(lapply(source_files, source))

args <- parse_args.cir()
n <- as.integer(get_arg.cir(args, "n", 500))
nsim <- as.integer(get_arg.cir(args, "nsim", 500))
example <- as.integer(get_arg.cir(args, "example", 1))
seed_fix <- as.integer(get_arg.cir(args, "seed", 3047))
slurm_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_)))
default_cores <- if (!is.na(slurm_cores) && slurm_cores > 0L) {
  slurm_cores
} else {
  max(1L, parallel::detectCores() - 1L)
}
cores <- as.integer(get_arg.cir(args, "cores", default_cores))
result_dir <- get_arg.cir(args, "result_dir", Sys.getenv("CIR_RESULT_DIR", unset = "Results"))
filename <- get_arg.cir(args, "filename",
                        sprintf("cir_example%d_n%d_nsim%d_seed%d", example, n, nsim, seed_fix))
cores <- max(1L, min(cores, nsim))

set.seed(seed_fix)
timepoints <- if (example == 1L) {
  1:10
} else if (example == 2L) {
  sort(stats::runif(10, 1, 10))
} else {
  stop("example must be 1 or 2.", call. = FALSE)
}

rho.tr   <- c(0.8, 0.25)
tau.tr   <- c(-1.8, 0.5)
gamma.tr <- c(0.1, -0.3)
eta.tr   <- c(-2.5, 0.3)

settings <- c(
  "all_correct",
  "sop_wrong_ps_correct_cen_correct",
  "sop_correct_ps_wrong_cen_wrong",
  "sop_wrong_ps_wrong_cen_wrong"
)

p <- length(rho.tr)
op.tr <- tau.tr
ps.tr <- gamma.tr
cen.tr <- eta.tr
names(op.tr) <- paste0("tau", seq_along(op.tr) - 1L)
names(ps.tr) <- paste0("gamma", seq_along(ps.tr) - 1L)
names(cen.tr) <- paste0("eta", seq_along(cen.tr) - 1L)
p_op <- length(op.tr)
p_ps <- length(ps.tr)
p_cen <- length(cen.tr)

point.est.dr <- array(NA_real_, c(nsim, p, length(settings)), dimnames = list(NULL, paste0("rho", 0:(p - 1)), settings))
sd.est.dr <- lcb.est.dr <- ucb.est.dr <- point.est.dr
point.est.mle <- array(NA_real_, c(nsim, p, length(settings)), dimnames = dimnames(point.est.dr))
sd.est.mle <- lcb.est.mle <- ucb.est.mle <- point.est.mle
point.est.op <- array(NA_real_, c(nsim, p_op, length(settings)), dimnames = list(NULL, names(op.tr), settings))
sd.est.op <- lcb.est.op <- ucb.est.op <- point.est.op
point.est.ps <- array(NA_real_, c(nsim, p_ps, length(settings)), dimnames = list(NULL, names(ps.tr), settings))
sd.est.ps <- lcb.est.ps <- ucb.est.ps <- point.est.ps
point.est.cen <- array(NA_real_, c(nsim, p_cen, length(settings)), dimnames = list(NULL, names(cen.tr), settings))
sd.est.cen <- lcb.est.cen <- ucb.est.cen <- point.est.cen

run_one_sim <- function(j) {
  dat <- Simul.Ex.cir(rho.tr, tau.tr, gamma.tr, timepoints, n, sample.seed = seed_fix+j)
  a <- dat$a
  x <- dat$x
  y <- dat$y
  delta_c <- dat$delta_c
  fun_tx <- dat$fun_tx
  x_wrong <- cbind(1, matrix(stats::runif(n * (ncol(x) - 1L), -2, 2), n, ncol(x) - 1L))
  x_wrong_ps <- cbind(1, matrix(stats::runif(n * (ncol(x) - 1L), -2, 2), n, ncol(x) - 1L))
  x_wrong_cen <- cbind(1, matrix(stats::runif(n * (ncol(x) - 1L), -2, 2), n, ncol(x) - 1L))

  pe_dr <- matrix(NA_real_, p, length(settings))
  sd_dr <- matrix(NA_real_, p, length(settings))
  lcb_dr <- matrix(NA_real_, p, length(settings))
  ucb_dr <- matrix(NA_real_, p, length(settings))
  pe_mle <- matrix(NA_real_, p, length(settings))
  sd_mle <- matrix(NA_real_, p, length(settings))
  lcb_mle <- matrix(NA_real_, p, length(settings))
  ucb_mle <- matrix(NA_real_, p, length(settings))
  pe_op <- matrix(NA_real_, p_op, length(settings))
  sd_op <- matrix(NA_real_, p_op, length(settings))
  lcb_op <- matrix(NA_real_, p_op, length(settings))
  ucb_op <- matrix(NA_real_, p_op, length(settings))
  pe_ps <- matrix(NA_real_, p_ps, length(settings))
  sd_ps <- matrix(NA_real_, p_ps, length(settings))
  lcb_ps <- matrix(NA_real_, p_ps, length(settings))
  ucb_ps <- matrix(NA_real_, p_ps, length(settings))
  pe_cen <- matrix(NA_real_, p_cen, length(settings))
  sd_cen <- matrix(NA_real_, p_cen, length(settings))
  lcb_cen <- matrix(NA_real_, p_cen, length(settings))
  ucb_cen <- matrix(NA_real_, p_cen, length(settings))
  failures <- character()

  for (s in seq_along(settings)) {
    setting <- settings[s]
    x_sop <- if (grepl("sop_wrong", setting)) x_wrong else x
    x_ps <- if (grepl("ps_wrong", setting)) x_wrong_ps else x
    x_cen <- if (grepl("cen_wrong", setting)) x_wrong_cen else x
    res <- tryCatch({

  ## 1. MLE
  mle_num <- MLEst.cir(
    a = a,
    fun_tx = fun_tx,
    x_sop = x_sop,
    y = y,
    delta_c = delta_c,
    rho.start = rep(0, p),
    tau.start = rep(0, length(tau.tr)),
    max.step = 1000,
    thres = 1e-6,
    timepoints = timepoints
  )

  ## 2. Pass the numerical MLE results to the analytic DR sandwich estimator
  DREst.cir.analytic(
    a = a,
    fun_tx = fun_tx,
    x_sop = x_sop,
    x_ps = x_ps,
    x_cen = x_cen,
    y = y,
    delta_c = delta_c,

    ## Key: directly use the numerical MLE fitted above
    mle.fit = mle_num,

    ## There is no need to provide tau.start to the analytic wrapper
    ## Initialize the DR estimator using rho.nuis from mle_num
    rho.start = NULL,

    ## The censoring model is still fitted by the analytic wrapper
    eta.start = rep(0, ncol(x_cen)),

    max.step = 1000,
    thres = 1e-6,
    sandwich = TRUE,
    timepoints = timepoints,

    dr.method = "checked",
    check_cdf = TRUE,
    cdf.tol = 0,
    score.tol = 1e-5
  )

}, error = function(e) {

  failures <<- c(
    failures,
    paste0(setting, ": ", conditionMessage(e))
  )

  NULL
})
    if (!is.null(res)) {
      pe_dr[, s] <- res$dr.point
      sd_dr[, s] <- res$sd.est
      lcb_dr[, s] <- res$conf.lower
      ucb_dr[, s] <- res$conf.upper
      pe_mle[, s] <- res$mle$point.est[seq_len(p)]
      sd_mle[, s] <- res$mle$se.est[seq_len(p)]
      lcb_mle[, s] <- res$mle$conf.lower[seq_len(p)]
      ucb_mle[, s] <- res$mle$conf.upper[seq_len(p)]
      op_idx <- p + seq_len(p_op)
      pe_op[, s] <- res$mle$point.est[op_idx]
      sd_op[, s] <- res$mle$se.est[op_idx]
      lcb_op[, s] <- res$mle$conf.lower[op_idx]
      ucb_op[, s] <- res$mle$conf.upper[op_idx]
      pe_ps[, s] <- res$ps$point.est
      sd_ps[, s] <- res$ps$se.est
      lcb_ps[, s] <- res$ps$conf.lower
      ucb_ps[, s] <- res$ps$conf.upper
      cen_se <- sqrt(diag(res$cen$cov))
      cen_z <- stats::qnorm(0.975)
      pe_cen[, s] <- res$cen$point.est
      sd_cen[, s] <- cen_se
      lcb_cen[, s] <- res$cen$point.est - cen_z * cen_se
      ucb_cen[, s] <- res$cen$point.est + cen_z * cen_se
    }
  }

  list(
    pe_dr = pe_dr,
    sd_dr = sd_dr,
    lcb_dr = lcb_dr,
    ucb_dr = ucb_dr,
    pe_mle = pe_mle,
    sd_mle = sd_mle,
    lcb_mle = lcb_mle,
    ucb_mle = ucb_mle,
    pe_op = pe_op,
    sd_op = sd_op,
    lcb_op = lcb_op,
    ucb_op = ucb_op,
    pe_ps = pe_ps,
    sd_ps = sd_ps,
    lcb_ps = lcb_ps,
    ucb_ps = ucb_ps,
    pe_cen = pe_cen,
    sd_cen = sd_cen,
    lcb_cen = lcb_cen,
    ucb_cen = ucb_cen,
    failures = failures
  )
}

message(sprintf("running %d simulations on %d core(s)", nsim, cores))
if (cores > 1L) {
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  work_dir <- getwd()
  parallel::clusterExport(cl, c("work_dir", "source_files"), envir = environment())
  parallel::clusterEvalQ(cl, {
    .libPaths(c("/home/yanjin41/R/4.3.1", .libPaths()))
    setwd(work_dir)
    invisible(lapply(source_files, source))
  })
  parallel::clusterExport(
    cl,
    c("n", "nsim", "seed_fix", "timepoints", "rho.tr", "tau.tr", "gamma.tr", "settings",
      "p", "p_op", "p_ps", "p_cen", "run_one_sim"),
    envir = environment()
  )
  res_list <- parallel::parLapply(cl, seq_len(nsim), run_one_sim)
} else {
  res_list <- lapply(seq_len(nsim), run_one_sim)
}

for (j in seq_len(nsim)) {
  point.est.dr[j, , ] <- res_list[[j]]$pe_dr
  sd.est.dr[j, , ] <- res_list[[j]]$sd_dr
  lcb.est.dr[j, , ] <- res_list[[j]]$lcb_dr
  ucb.est.dr[j, , ] <- res_list[[j]]$ucb_dr
  point.est.mle[j, , ] <- res_list[[j]]$pe_mle
  sd.est.mle[j, , ] <- res_list[[j]]$sd_mle
  lcb.est.mle[j, , ] <- res_list[[j]]$lcb_mle
  ucb.est.mle[j, , ] <- res_list[[j]]$ucb_mle
  point.est.op[j, , ] <- res_list[[j]]$pe_op
  sd.est.op[j, , ] <- res_list[[j]]$sd_op
  lcb.est.op[j, , ] <- res_list[[j]]$lcb_op
  ucb.est.op[j, , ] <- res_list[[j]]$ucb_op
  point.est.ps[j, , ] <- res_list[[j]]$pe_ps
  sd.est.ps[j, , ] <- res_list[[j]]$sd_ps
  lcb.est.ps[j, , ] <- res_list[[j]]$lcb_ps
  ucb.est.ps[j, , ] <- res_list[[j]]$ucb_ps
  point.est.cen[j, , ] <- res_list[[j]]$pe_cen
  sd.est.cen[j, , ] <- res_list[[j]]$sd_cen
  lcb.est.cen[j, , ] <- res_list[[j]]$lcb_cen
  ucb.est.cen[j, , ] <- res_list[[j]]$ucb_cen
  if (length(res_list[[j]]$failures) > 0L) {
    message("simulation ", j, " failures: ", paste(res_list[[j]]$failures, collapse = " | "))
  }
}

summarize_est <- function(point, se, lcb, ucb, truth) {
  S <- dim(point)[3]
  p <- dim(point)[2]
  slice_setting <- function(arr, s) {
    matrix(arr[, , s], ncol = p, dimnames = list(NULL, dimnames(point)[[2]]))
  }
  bias_mean <- bias_se <- ratio_se <- coverage <- matrix(NA_real_, S, p, dimnames = list(dimnames(point)[[3]], dimnames(point)[[2]]))
  for (s in seq_len(S)) {
    pe_all <- slice_setting(point, s)
    se_all <- slice_setting(se, s)
    lo_all <- slice_setting(lcb, s)
    hi_all <- slice_setting(ucb, s)
    keep <- complete.cases(pe_all, se_all, lo_all, hi_all)
    pe <- pe_all[keep, , drop = FALSE]
    see <- se_all[keep, , drop = FALSE]
    lo <- lo_all[keep, , drop = FALSE]
    hi <- hi_all[keep, , drop = FALSE]
    if (nrow(pe) == 0L) next
    bias <- sweep(pe, 2, truth, "-")
    bias_mean[s, ] <- colMeans(bias)
    bias_se[s, ] <- apply(bias, 2, stats::sd) / sqrt(nrow(pe))
    mc_se <- apply(pe, 2, stats::sd)
    ratio_se[s, ] <- colMeans(see) / mc_se
    coverage[s, ] <- sapply(seq_len(p), function(k) mean(lo[, k] < truth[k] & hi[, k] > truth[k]))
  }
  list(bias_mean = bias_mean, bias_se = bias_se, ratio_se = ratio_se, coverage = coverage)
}

summaries <- list(
  MLE = summarize_est(point.est.mle, sd.est.mle, lcb.est.mle, ucb.est.mle, rho.tr),
  DR = summarize_est(point.est.dr, sd.est.dr, lcb.est.dr, ucb.est.dr, rho.tr),
  OP = summarize_est(point.est.op, sd.est.op, lcb.est.op, ucb.est.op, op.tr),
  PS = summarize_est(point.est.ps, sd.est.ps, lcb.est.ps, ucb.est.ps, ps.tr),
  CEN = summarize_est(point.est.cen, sd.est.cen, lcb.est.cen, ucb.est.cen, cen.tr)
)

array_to_long.cir <- function(point, se, lcb, ucb, truth, method) {
  nsim_local <- dim(point)[1]
  p_local <- dim(point)[2]
  S_local <- dim(point)[3]
  param_names <- dimnames(point)[[2]]
  setting_names <- dimnames(point)[[3]]
  out <- vector("list", nsim_local * p_local * S_local)
  idx <- 0L
  for (s in seq_len(S_local)) {
    for (j in seq_len(nsim_local)) {
      for (k in seq_len(p_local)) {
        idx <- idx + 1L
        out[[idx]] <- data.frame(
          method = method,
          sim = j,
          setting = setting_names[s],
          parameter = param_names[k],
          truth = truth[k],
          estimate = point[j, k, s],
          bias = point[j, k, s] - truth[k],
          se = se[j, k, s],
          conf.lower = lcb[j, k, s],
          conf.upper = ucb[j, k, s],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, out)
}

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
path.save <- file.path(result_dir, filename)
if (!dir.exists(path.save)) dir.create(path.save, recursive = TRUE)
raw_results <- rbind(
  array_to_long.cir(point.est.mle, sd.est.mle, lcb.est.mle, ucb.est.mle, rho.tr, "MLE"),
  array_to_long.cir(point.est.dr, sd.est.dr, lcb.est.dr, ucb.est.dr, rho.tr, "DR"),
  array_to_long.cir(point.est.op, sd.est.op, lcb.est.op, ucb.est.op, op.tr, "OP"),
  array_to_long.cir(point.est.ps, sd.est.ps, lcb.est.ps, ucb.est.ps, ps.tr, "PS"),
  array_to_long.cir(point.est.cen, sd.est.cen, lcb.est.cen, ucb.est.cen, cen.tr, "CEN")
)
write.csv(raw_results, file = file.path(path.save, "raw_results.csv"), row.names = FALSE)
save.image(file = file.path(path.save, "Simu.RData"))
sink(file.path(path.save, "simu.out"))
date()
cat("code_dir\n")
print(code_dir)
cat("result_dir\n")
print(result_dir)
cat("n\n")
print(n)
cat("nsim\n")
print(nsim)
cat("cores\n")
print(cores)
cat("example\n")
print(example)
cat("timepoints\n")
print(timepoints)
cat("settings\n")
print(settings)
cat("summaries\n")
print(summaries)
sink()
