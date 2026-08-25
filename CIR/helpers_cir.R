################################################################################
############################ CIR chapter 4 helpers #############################
################################################################################

expit.cir <- function(logodds) {
  1 / (1 + exp(-logodds))
}

as_matrix.cir <- function(x, n = NULL, name = "x") {
  if (is.null(dim(x))) {
    out <- matrix(x, ncol = 1)
  } else {
    out <- as.matrix(x)
  }
  if (!is.null(n) && nrow(out) != n) {
    stop(name, " must have ", n, " rows.", call. = FALSE)
  }
  out
}

param_dim.cir <- function(x) {
  if (is.null(dim(x))) 1 else ncol(x)
}

fun_tx_dim.cir <- function(fun_tx) {
  param_dim.cir(fun_tx[[1]])
}

check_fun_tx.cir <- function(fun_tx, n = NULL, m = NULL) {
  if (!is.list(fun_tx)) {
    stop("fun_tx must be a list.", call. = FALSE)
  }
  if (!is.null(m) && length(fun_tx) != m + 1) {
    stop("fun_tx must have length m + 1.", call. = FALSE)
  }
  if (length(fun_tx) < 2) {
    stop("fun_tx must contain at least two time design matrices.", call. = FALSE)
  }
  p <- fun_tx_dim.cir(fun_tx)
  for (j in seq_along(fun_tx)) {
    xj <- as_matrix.cir(fun_tx[[j]], name = paste0("fun_tx[[", j, "]]"))
    if (!is.null(n) && nrow(xj) != n) {
      stop("fun_tx[[", j, "]] must have ", n, " rows.", call. = FALSE)
    }
    if (ncol(xj) != p) {
      stop("All fun_tx elements must have the same number of columns.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

resolve_timepoints.cir <- function(y, fun_tx = NULL, timepoints = NULL) {
  if (is.null(timepoints) && !is.null(fun_tx)) {
    timepoints <- attr(fun_tx, "timepoints", exact = TRUE)
  }
  if (is.null(timepoints)) {
    timepoints <- sort(unique(y))
  } else {
    timepoints <- sort(unique(timepoints))
  }
  if (any(!y %in% timepoints)) {
    stop("All observed y values must be included in timepoints.", call. = FALSE)
  }
  if (!is.null(fun_tx) && length(fun_tx) != length(timepoints) + 1L) {
    stop("fun_tx length must equal length(timepoints) + 1.", call. = FALSE)
  }
  timepoints
}

fun_tx_matrix.cir <- function(fun_tx, index, n = NULL) {
  as_matrix.cir(fun_tx[[index]], n = n, name = paste0("fun_tx[[", index, "]]"))
}

fun_tx_row_stack.cir <- function(fun_tx, i, m) {
  do.call(rbind, lapply(seq_len(m), function(k) {
    as.numeric(fun_tx[[k]][i, , drop = TRUE])
  }))
}

safe_solve.cir <- function(mat, rhs = NULL, ridge = 1e-8) {
  mat <- as.matrix(mat)
  ans <- tryCatch({
    if (is.null(rhs)) solve(mat) else solve(mat, rhs)
  }, error = function(e) NULL)
  if (!is.null(ans)) {
    return(ans)
  }
  mat2 <- mat + diag(ridge, nrow(mat))
  if (is.null(rhs)) solve(mat2) else solve(mat2, rhs)
}

parse_args.cir <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0) {
    return(list())
  }
  kv <- strsplit(sub("^--", "", args), "=", fixed = TRUE)
  setNames(lapply(kv, `[`, 2L), vapply(kv, `[`, "", 1L))
}

get_arg.cir <- function(arg_list, name, default = NULL) {
  val <- arg_list[[name]]
  if (is.null(val) || is.na(val) || identical(val, "")) default else val
}
