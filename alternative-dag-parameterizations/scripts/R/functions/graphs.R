# Graph definitions use d-separation, rather than sampled coefficient support.
imap_closure <- function(A) {
  R <- A != 0
  for (k in seq_len(nrow(A))) R <- R | outer(R[, k], R[k, ], "&")
  R
}

# Ancestral moralization is independent of the Gaussian coefficient checks.
imap_dseparated <- function(A, u, v, C = integer()) {
  A <- A != 0
  R <- imap_closure(A); endpoints <- c(u, v, C)
  anc <- sort(unique(c(endpoints, which(rowSums(R[, endpoints, drop = FALSE]) > 0))))
  moral <- A | t(A)
  for (j in anc) {
    P <- which(A[, j]); moral[P, P] <- TRUE
  }
  allowed <- setdiff(anc, C); visited <- u; frontier <- u
  while (length(frontier)) {
    next_nodes <- setdiff(allowed[colSums(moral[frontier, allowed, drop = FALSE]) > 0], visited)
    visited <- union(visited, next_nodes); frontier <- next_nodes
  }
  !(v %in% visited)
}

imap_ordered_graph <- function(A, order) {
  G <- matrix(FALSE, nrow(A), ncol(A))
  for (k in seq_along(order)) {
    j <- order[k]; P <- order[seq_len(k - 1L)]
    for (u in P) G[u, j] <- !imap_dseparated(A, u, j, setdiff(P, u))
  }
  G
}

imap_topological_order <- function(A) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A), all(A %in% c(0, 1)))
  remaining <- seq_len(nrow(A)); order <- integer()
  while (length(remaining)) {
    roots <- remaining[colSums(A[remaining, remaining, drop = FALSE]) == 0]
    if (!length(roots)) stop("The adjacency matrix contains a directed cycle")
    # A deterministic tie-break keeps results reproducible.
    j <- min(roots); order <- c(order, j); remaining <- setdiff(remaining, j)
  }
  order
}

imap_orders <- function(A, d, y) {
  topo <- imap_topological_order(A)
  R <- imap_closure(A)
  stopifnot(R[d, y])
  M <- topo[topo %in% which(R[d, ] & R[, y])]; S <- c(d, M, y)
  ancestors <- setdiff(which(rowSums(R[, S, drop = FALSE]) > 0), S)
  delta1 <- topo[topo %in% ancestors]
  delta3 <- setdiff(topo, c(delta1, S))
  original <- c(delta1, d, M, y, delta3)
  optimal <- c(delta1, d, y, M, delta3)
  O <- setdiff(which(rowSums(A[, c(M, y), drop = FALSE]) > 0), c(d, M))
  list(original = as.integer(original), optimal = as.integer(optimal), M = M, O = O)
}
