# Independent finite graphical and covariance benchmarks.
imap_validate <- function(seed = 20261010L, tolerance = 1e-9) {
  set.seed(seed)
  errors <- c(covariance = 0, regression = 0, round_trip = 0, effect = 0,
              diamond = 0, constraints = 0)
  q <- 5L; possible <- which(upper.tri(matrix(FALSE, q, q)), arr.ind = TRUE)
  pair_count <- 0L; swaps <- 0L; graph_count <- 0L
  for (mask in 0:(2^nrow(possible) - 1L)) {
    A <- matrix(FALSE, q, q)
    bits <- as.logical(intToBits(mask)[seq_len(nrow(possible))]); A[possible] <- bits
    R <- imap_closure(A); pairs <- which(R, arr.ind = TRUE)
    graph_count <- graph_count + 1L
    L <- diag(q); L[A] <- runif(sum(A), -0.8, 0.8)
    v <- runif(q, 0.5, 2); Sigma <- imap_covariance(L, v)
    for (k in seq_len(nrow(pairs))) {
      d <- pairs[k, 1]; y <- pairs[k, 2]; o <- imap_orders(A, d, y)
      G <- imap_ordered_graph(A, o$optimal)
      if (!setequal(which(G[, y]), c(d, o$O)))
        stop(paste('Parent mismatch: mask', mask, 'd', d, 'y', y,
                   'got', paste(which(G[,y]),collapse=','),
                   'expected', paste(sort(c(d,o$O)),collapse=',')))
      p <- list(L = L, v = v, ordering = o$original)
      for (m in rev(o$M)) p <- imap_transpose(p, m, y)
      stopifnot(identical(p$ordering, o$optimal), all(p$v > 0))
      reference <- imap_cholesky(Sigma, o$optimal)
      errors['covariance'] <- max(errors['covariance'], abs(imap_covariance(p$L, p$v) - Sigma))
      errors['regression'] <- max(errors['regression'], abs(p$L - reference$L), abs(p$v - reference$v))
      # The graph is defined by d-separation; accidental sampled cancellation is not a graph definition.
      stopifnot(all(abs(p$L[!G & row(p$L) != col(p$L)]) < tolerance))
      total <- solve(L)[d, y]
      errors['effect'] <- max(errors['effect'], abs(-p$L[d, y] - total))
      back <- p
      for (m in o$M) back <- imap_transpose(back, y, m)
      errors['round_trip'] <- max(errors['round_trip'], abs(back$L - L), abs(back$v - v))
      pair_count <- pair_count + 1L; swaps <- swaps + length(o$M)
    }
  }
  # Four-node diamond, including non-unit error variances.
  a <- -0.7; b <- 0.4; c <- -0.6; e <- 0.8
  v <- c(1.2, 0.7, 1.5, 0.9); L <- diag(4)
  L[1, 2] <- a; L[1, 3] <- b; L[2, 4] <- c; L[3, 4] <- e
  p <- list(L = L, v = v, ordering = 1:4)
  p <- imap_transpose(imap_transpose(p, 3, 4), 2, 4)
  s3 <- v[4] + e^2 * v[3]; s <- s3 + c^2 * v[2]; t0 <- a*c + b*e
  expected <- rbind(c(1,0,0,0), c(a-t0*c*v[2]/s,1,0,c*v[2]/s),
                    c(b-b*e^2*v[3]/s3,c*e*v[3]/s3,1,e*v[3]/s3),
                    c(-t0,0,0,1))
  errors['diamond'] <- max(abs(t(p$L) - expected),
                          abs(p$v - c(v[1], v[2]*s3/s, v[3]*v[4]/s3, s)))
  # General-variance implicit constraints; never divide by a coefficient.
  z <- p$L; w <- p$v
  f <- z[4,2]; g <- z[4,3]
  constraints <- c(z[2,3]*(w[2]+f^2*w[4])-g*f*w[4],
                   w[3]*(z[1,4]*w[2]+f*w[4]*z[1,2])+g*w[4]*w[2]*z[1,3])
  errors['constraints'] <- max(abs(constraints))
  # Reverse an unconstrained I-MAP: the same equations detect nonmembers.
  for (k in 1:100) {
    z <- diag(4); z[1,2] <- runif(1,-1,1); z[1,3] <- runif(1,-1,1)
    z[1,4] <- runif(1,-1,1); z[2,3] <- runif(1,-1,1)
    z[4,2] <- runif(1,-1,1); z[4,3] <- runif(1,-1,1); w <- runif(4,0.5,2)
    f <- z[4,2]; g <- z[4,3]; s2 <- w[2]+f^2*w[4]; t4 <- w[4]*w[2]/s2
    A14 <- (z[1,4]*w[2]+f*w[4]*z[1,2])/s2
    back <- imap_transpose(imap_transpose(list(L=z,v=w,ordering=c(1L,4L,2L,3L)),4,2),4,3)
    errors['constraints'] <- max(errors['constraints'],
      abs(back$L[2,3]-(z[2,3]-g*f*w[4]/s2)),
      abs(back$L[1,4]-(w[3]*A14+g*t4*z[1,3])/(w[3]+g^2*t4)))
  }
  # Cancellation: graph has causal paths, but their sum is exactly zero.
  cancel <- diag(4); cancel[1,2] <- -1; cancel[1,3] <- -1
  cancel[2,4] <- -1; cancel[3,4] <- 1
  pc <- imap_transpose(imap_transpose(list(L=cancel,v=rep(1,4),ordering=1:4),3,4),2,4)
  stopifnot(pc$L[1,4] == 0)
  bad <- tryCatch(imap_transpose(list(L=L,v=v,ordering=1:4),1,4),error=function(e) e)
  stopifnot(inherits(bad, 'error'), all(errors < tolerance))
  list(summary = data.frame(check = c('DAGs in fixed topological order', 'Ancestor pairs',
         'Forward adjacent swaps', 'Unconstrained reverse examples'),
         cases = c(graph_count,pair_count,swaps,100)),
       errors = data.frame(check = names(errors), maximum_absolute_error = unname(errors)),
       diamond = data.frame(parameter = c('total effect','outcome residual variance'),
                            value = c(-p$L[1,4],p$v[4])), seed = seed, tolerance = tolerance)
}
