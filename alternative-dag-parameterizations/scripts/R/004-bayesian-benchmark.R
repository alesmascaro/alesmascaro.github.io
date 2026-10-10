# Q004: exact three-graph baseline/normalized-moment evidence and posterior sampling.
# Known zero means; synthetic chain with B12=.7, B23=.9, variances 1, n=120.
# Compatible baseline a=8, U=I; 8,000 conditional posterior draws, seed 20261012.
source("scripts/R/functions/load.R")
bayesian_results <- imap_validate_bayesian()
if (sys.nframe() == 0L) {
  dir.create("results/tables",recursive=TRUE,showWarnings=FALSE)
  write.csv(bayesian_results$scores,"results/tables/004-bayesian-scores.csv",row.names=FALSE)
  write.csv(bayesian_results$monte_carlo,"results/tables/004-bayesian-moments.csv",row.names=FALSE)
  print(bayesian_results$scores); print(bayesian_results$moments)
  print(bayesian_results$monte_carlo); print(bayesian_results$diagnostics)
  writeLines(c(paste("seed:",bayesian_results$seed),paste("samples:",bayesian_results$samples),
    paste("n:",bayesian_results$n),paste("tolerance:",bayesian_results$tolerance),
    capture.output(sessionInfo())),"results/004-bayesian-session.txt")
}
