# Run from the repository root. This loader defines functions without running analyses.
for (module in c("graphs.R", "parameters.R", "frugal.R", "dag-wishart.R",
                 "adjustment.R", "plots.R", "validate-imap.R", "validate-workflows.R",
                 "effect-testing.R")) {
  source(file.path("scripts/R/functions", module), local = TRUE)
}
rm(module)
