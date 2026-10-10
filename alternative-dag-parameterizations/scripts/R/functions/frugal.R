# Explicit Gaussian marginal-and-association coordinates for the three-node chain.
imap_chain_frugal <- function(b, c, v) {
  stopifnot(length(v) == 3L, all(is.finite(v)), all(v > 0))
  s <- c^2 * v[2] + v[3]
  list(v1 = v[1], b = b, v2 = v[2], tau = b * c, s = s,
       rho = c * sqrt(v[2] / s))
}

# The six-coordinate ambient model permits a direct 1 -> 3 edge.
imap_frugal_parameters <- function(z) {
  stopifnot(all(c(z$v1, z$v2, z$s) > 0), abs(z$rho) < 1)
  alpha <- z$rho * sqrt(z$s / z$v2)
  L <- diag(3)
  L[1, 2] <- -z$b; L[2, 3] <- -alpha
  L[1, 3] <- -(z$tau - alpha * z$b)
  list(L = L, v = c(z$v1, z$v2, z$s * (1 - z$rho^2)), ordering = 1:3)
}

imap_frugal_chain_inverse <- function(z, tolerance = 1e-10) {
  p <- imap_frugal_parameters(z)
  if (abs(p$L[1, 3]) > tolerance) stop("Frugal coordinates violate the chain constraint")
  list(b = -p$L[1, 2], c = -p$L[2, 3], v = p$v)
}
