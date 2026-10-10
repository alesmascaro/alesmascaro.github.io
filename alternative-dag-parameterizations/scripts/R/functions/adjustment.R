# Sufficient back-door test, complete within the stipulated pretreatment candidates.
# This function does not implement the generalized adjustment criterion for all sets.
imap_valid_pretreatment_adjustment <- function(A, d, y, Z) {
  imap_topological_order(A)
  R <- imap_closure(A)
  if (any(Z %in% c(d, y, which(R[d, ])))) return(FALSE)
  backdoor <- A; backdoor[d, ] <- FALSE
  imap_dseparated(backdoor, d, y, Z)
}

# The four synthetic graphs in the archived two-step proof-of-concept script.
imap_adjustment_example_graphs <- function() {
  A <- matrix(FALSE, 6, 6)
  A[rbind(c(1,2), c(1,4), c(4,5), c(5,6), c(4,6), c(2,5))] <- TRUE
  graphs <- rep(list(A), 4)
  graphs[[1]][3,6] <- TRUE; graphs[[2]][6,3] <- TRUE
  graphs[[4]][2,5] <- FALSE; graphs[[4]][5,2] <- TRUE
  names(graphs) <- paste0("D", 1:4)
  graphs
}
