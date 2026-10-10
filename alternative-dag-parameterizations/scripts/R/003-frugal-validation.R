# Q005-Q006: sparse frugal coordinates, arbitrary node labels and adjustment validity.
# All computations use population parameters. Run from the repository root.
source("scripts/R/functions/load.R")
frugal_results <- imap_validate_frugal()
if (sys.nframe() == 0L) {
  dir.create("results/tables",recursive=TRUE,showWarnings=FALSE)
  write.csv(frugal_results$errors,"results/tables/003-frugal-errors.csv",row.names=FALSE)
  write.csv(frugal_results$validity,"results/tables/003-adjustment-validity.csv",row.names=FALSE)
  print(frugal_results$summary); print(frugal_results$errors); print(frugal_results$validity)
  writeLines(c(paste("seed:",frugal_results$seed),paste("tolerance:",frugal_results$tolerance),
    capture.output(sessionInfo())),"results/003-frugal-session.txt")
}
