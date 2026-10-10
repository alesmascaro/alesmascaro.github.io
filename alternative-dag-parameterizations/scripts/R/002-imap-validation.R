# Q001-Q003: main validation script. Run from the repository root.
source("scripts/R/functions/load.R")
imap_results <- imap_validate()

if (sys.nframe() == 0L) {
  dir.create('results/tables', recursive = TRUE, showWarnings = FALSE)
  write.csv(imap_results$summary, 'results/tables/002-imap-cases.csv', row.names = FALSE)
  write.csv(imap_results$errors, 'results/tables/002-imap-errors.csv', row.names = FALSE)
  print(imap_results$summary); print(imap_results$errors)
  writeLines(c(paste('seed:', imap_results$seed), paste('tolerance:', imap_results$tolerance),
               capture.output(sessionInfo())), 'results/002-imap-session.txt')
}
