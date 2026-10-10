# Q004/Q006: exact 25-DAG comparison of causal-effect hypothesis probabilities.
# Synthetic known-zero-mean chain/collider, B=.8, errors 1; equal graph prior mass.
# Baseline a=8, U=I3; active pi_M=tau^2*pi_0/E_0[tau^2], unchanged null priors.
# 200 replicates, nested n=25,50,100,250,500,1000; seed 20261013.
source("scripts/R/functions/load.R")
effect_moment_checks <- imap_validate_effect_moments()
effect_testing_results <- imap_simulate_effect_testing()
if (sys.nframe() == 0L) {
  dir.create("results/tables", recursive=TRUE, showWarnings=FALSE)
  dir.create("results/figures", recursive=TRUE, showWarnings=FALSE)
  write.csv(effect_testing_results$raw, "results/tables/005-effect-testing-replicates.csv", row.names=FALSE)
  write.csv(effect_testing_results$summary, "results/tables/005-effect-testing-summary.csv", row.names=FALSE)
  write.csv(effect_moment_checks$checks, "results/tables/005-effect-moment-checks.csv", row.names=FALSE)
  png("results/figures/005-effect-testing.png", width=1600, height=1150, res=180)
  imap_plot_effect_testing(effect_testing_results); dev.off()
  pdf("results/figures/005-effect-testing.pdf", width=9, height=6.5)
  imap_plot_effect_testing(effect_testing_results); dev.off()
  print(effect_testing_results$summary)
  writeLines(c(paste("seed:",effect_testing_results$seed),
    paste("replicates:",effect_testing_results$replicates),
    paste("sample sizes:",paste(effect_testing_results$sizes,collapse=",")),
    "a=8; U=I3; all 25 DAGs; uniform graph mass; B=.8; error variances=1; known zero means",
    capture.output(sessionInfo())), "results/005-effect-testing-session.txt")
}
