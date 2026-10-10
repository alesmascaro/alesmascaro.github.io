# Exact baseline in the book's L = I - B convention; known zero means.
dagwishart_shapes <- function(A, a) {
  q <- nrow(A)
  stopifnot(a > q - 1)
  a - q + 2 * colSums(A) + 3
}

dagwishart_blocks <- function(A, alpha, U) {
  imap_topological_order(A)
  q <- nrow(A)
  stopifnot(length(alpha) == q, all(dim(U) == c(q, q)), isSymmetric(U))
  chol(U) # Verifies positive definiteness.
  lapply(seq_len(q), function(j) {
    P <- which(A[, j] != 0); p <- length(P)
    r <- (alpha[j] - p) / 2 - 1
    if (r <= 0) stop("DAG-Wishart shape is not proper")
    if (p) {
      K <- U[P, P, drop = FALSE]; Kinv <- chol2inv(chol(K))
      u <- U[P, j]; mu <- -as.vector(Kinv %*% u)
      s <- U[j, j] - sum(u * (Kinv %*% u))
      logdet <- as.numeric(determinant(K, logarithm = TRUE)$modulus)
    } else {
      Kinv <- matrix(numeric(), 0, 0); mu <- numeric(); s <- U[j, j]; logdet <- 0
    }
    list(P = P, p = p, r = r, rate = s / 2, mu = mu,
         Kinv = Kinv, logdet = logdet)
  })
}

dagwishart_log_normalizer <- function(A, alpha, U) {
  sum(vapply(dagwishart_blocks(A, alpha, U), function(z)
    z$p / 2 * log(2 * pi) - z$logdet / 2 + lgamma(z$r) - z$r * log(z$rate), 0.0))
}

dagwishart_log_marginal <- function(X, A, a, U) {
  n <- nrow(X); q <- ncol(X); alpha <- dagwishart_shapes(A, a)
  -n * q / 2 * log(2 * pi) +
    dagwishart_log_normalizer(A, alpha + n, U + crossprod(X)) -
    dagwishart_log_normalizer(A, alpha, U)
}

dagwishart_sample <- function(S, A, alpha, U) {
  stopifnot(S >= 1, S == as.integer(S))
  q <- nrow(A); L <- array(0, c(q, q, S)); v <- matrix(0, S, q)
  blocks <- dagwishart_blocks(A, alpha, U)
  for (j in seq_len(q)) {
    z <- blocks[[j]]; L[j, j, ] <- 1
    v[, j] <- 1 / rgamma(S, shape = z$r, rate = z$rate)
    if (z$p) {
      noise <- matrix(rnorm(S * z$p), S, z$p) %*% chol(z$Kinv)
      coefficients <- sweep(noise * sqrt(v[, j]), 2, z$mu, "+")
      for (k in seq_along(z$P)) L[z$P[k], j, ] <- coefficients[, k]
    }
  }
  list(L = L, v = v)
}

# Moments of one structural coefficient B[u,j] under a node's N/IG block.
dagwishart_coefficient_moments <- function(block, u) {
  k <- match(u, block$P)
  stopifnot(!is.na(k), block$r > 2)
  mu <- -block$mu[k]; variance <- block$rate / (block$r - 1) * block$Kinv[k, k]
  fourth_central <- 3 * block$rate^2 / ((block$r - 1) * (block$r - 2)) * block$Kinv[k, k]^2
  c(first = mu, second = mu^2 + variance, third = mu^3 + 3 * mu * variance,
    fourth = mu^4 + 6 * mu^2 * variance + fourth_central)
}

imap_normalized_weights <- function(logw) {
  if (!any(is.finite(logw))) stop("No finite importance weights")
  w <- exp(logw - max(logw)); w / sum(w)
}
