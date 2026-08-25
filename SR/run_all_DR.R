rm(list = ls())
# run_all.R
date_tag <- format(Sys.Date(), "%m%d")
# seed_fix <- 3047
seed_fix <- 20250820

make_fname <- function(example, n) {
  sprintf("example%d_n%d_seed%d_%s", example, n, seed_fix, date_tag)
}

combos <- expand.grid(example = c(1L, 2L), n = c(500L, 1000L), KEEP.OUT.ATTRS = FALSE)

for (i in seq_len(nrow(combos))) {
  ex <- combos$example[i]
  nn <- combos$n[i]
  fname <- make_fname(ex, nn)
  
  cmd_args <- c("run_simul_DR.R",
                sprintf("--n=%d", nn),
                sprintf("--example=%d", ex),
                sprintf("--filename=%s", fname))
  
  message("--------------------------------------------------")
  message(sprintf("Calling: Rscript %s", paste(cmd_args, collapse = " ")))
  status <- system2("Rscript", args = cmd_args)
  if (status != 0) {
    warning(sprintf("Run failed for example=%d, n=%d (status %d)", ex, nn, status))
  }
}
message("All runs attempted.")

