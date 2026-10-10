# Workflow demonstration only; not an I-MAP result.
# Run from the repository root with Rscript scripts/R/001-monte-carlo.R.
# Source this same file from Quarto to compute the book's results.
simulate_means <- function(sample_sizes = c(10L, 100L, 1000L),
                           repetitions = 2000L, seed = 20261010L) {
  stopifnot(all(sample_sizes > 0), repetitions > 1)
  set.seed(seed)
  do.call(rbind, lapply(sample_sizes, function(n) {
    means <- replicate(repetitions, mean(rnorm(n)))
    data.frame(sample_size = n, repetitions = repetitions,
               empirical_mean = mean(means),
               empirical_sd = sd(means), theoretical_sd = 1 / sqrt(n))
  }))
}

mc_results <- simulate_means()

if (sys.nframe() == 0L) {
  dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
  write.csv(mc_results, "results/tables/001-monte-carlo.csv", row.names = FALSE)
  print(mc_results)
  writeLines(capture.output(sessionInfo()), "results/001-r-session.txt")
}
