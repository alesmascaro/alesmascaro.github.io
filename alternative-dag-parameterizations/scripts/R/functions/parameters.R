# Gaussian factorization utilities. L = I - B; A[u,v] denotes u -> v.
imap_covariance <- function(L, v) {
  K <- solve(L)
  t(K) %*% diag(v, nrow = length(v)) %*% K
}

imap_cholesky <- function(Sigma, ordering) {
  q <- nrow(Sigma); L <- diag(q); v <- numeric(q)
  for (k in seq_along(ordering)) {
    j <- ordering[k]; P <- ordering[seq_len(k - 1L)]
    if (length(P)) {
      beta <- solve(Sigma[P, P, drop = FALSE], Sigma[P, j])
      L[P, j] <- -beta
      v[j] <- Sigma[j, j] - sum(Sigma[j, P] * beta)
    } else v[j] <- Sigma[j, j]
  }
  list(L = L, v = v, ordering = ordering)
}

imap_transpose <- function(parameters, i, j) {
  L <- parameters$L; v <- parameters$v; order <- parameters$ordering
  ki <- match(i, order); kj <- match(j, order)
  if (is.na(ki) || is.na(kj) || kj != ki + 1L)
    stop("Transposition requires adjacent nodes i,j in the current ordering")
  a <- L[i, j]; vi <- v[i]; vj <- v[j]; s <- vj + a^2 * vi
  L[, j] <- L[, j] - a * L[, i]
  L[, i] <- L[, i] + (a * vi / s) * L[, j]
  v[i] <- vi * vj / s; v[j] <- s # Both use saved OLD variances.
  order[c(ki, kj)] <- order[c(kj, ki)]
  list(L = L, v = v, ordering = order)
}

# Reviewed counterpart of the archived DAGtoIMAP: explicitly checks the ordering.
imap_transform <- function(L, v, A, d, y) {
  stopifnot(nrow(L) == nrow(A), length(v) == nrow(A), all(v > 0),
            all(diag(L) == 1), all(L[!as.logical(A) & row(L) != col(L)] == 0))
  orders <- imap_orders(A, d, y)
  out <- list(L = L, v = v, ordering = orders$original)
  for (m in rev(orders$M)) out <- imap_transpose(out, m, y)
  out
}

# Reviewed counterpart of DAGtoMCM; replaces treatment by a constant in the SEM.
imap_intervention <- function(L, v, d, y) {
  stopifnot(d != y, all(v > 0))
  cut <- L; cut[, d] <- 0; cut[d, d] <- 1
  K <- solve(cut)
  list(effect = K[d, y], variance = sum(v[-d] * K[-d, y]^2))
}
