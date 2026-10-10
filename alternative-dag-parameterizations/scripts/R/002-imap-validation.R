# Q001-Q003: algebraic and finite graphical checks for the corrected I-MAP notes.
# Base R only. Run from the project root; Quarto sources the same file.
# Adjacency matrices use A[u,v] for u -> v; L = I - B.
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

imap_closure <- function(A) {
  R <- A
  for (k in seq_len(nrow(A))) R <- R | outer(R[, k], R[k, ], "&")
  R
}

# Ancestral moralization is independent of the Gaussian coefficient checks.
imap_dseparated <- function(A, u, v, C = integer()) {
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

imap_orders <- function(A, d, y) {
  R <- imap_closure(A)
  stopifnot(R[d, y])
  M <- which(R[d, ] & R[, y]); S <- c(d, M, y)
  delta1 <- setdiff(which(rowSums(R[, S, drop = FALSE]) > 0), S)
  # Enumeration inputs have natural topological order 1,...,q.
  original <- c(delta1, d, M, y, setdiff(seq_len(nrow(A)), c(delta1, S)))
  optimal <- c(delta1, d, y, M, setdiff(seq_len(nrow(A)), c(delta1, S)))
  O <- setdiff(which(rowSums(A[, c(M, y), drop = FALSE]) > 0), c(d, M))
  list(original = as.integer(original), optimal = as.integer(optimal), M = M, O = O)
}

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
  # Corrected four-node diamond, including non-unit error variances.
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

imap_results <- imap_validate()
imap_plot_diamond <- function() {
  old <- par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 0.5))
  on.exit(par(old))
  positions <- list(rbind(c(.5,.9),c(.15,.5),c(.85,.5),c(.5,.1)),
                    rbind(c(.1,.85),c(.45,.5),c(.45,.1),c(.85,.85)))
  edges <- list(rbind(c(1,2),c(1,3),c(2,4),c(3,4)),
                rbind(c(1,4),c(1,2),c(4,2),c(1,3),c(4,3),c(2,3)))
  titles <- c('Original causal DAG', 'Ordered I-MAP')
  for (k in 1:2) {
    xy <- positions[[k]]; plot.new(); plot.window(c(0,1), c(0,1), asp = 1)
    title(titles[k], cex.main = 1)
    for (r in seq_len(nrow(edges[[k]]))) {
      u <- edges[[k]][r,1]; v <- edges[[k]][r,2]
      direction <- xy[v,]-xy[u,]; direction <- direction/sqrt(sum(direction^2))
      start <- xy[u,]+.065*direction; end <- xy[v,]-.065*direction
      color <- if (k == 2 && u == 1 && v == 4) 'steelblue4' else 'gray35'
      arrows(start[1],start[2],end[1],end[2],length=.09,col=color,lwd=1.5)
    }
    points(xy,pch=21,bg='white',cex=3.5,lwd=1.3)
    text(xy,labels=1:4,cex=.9)
    if (k == 2) text(.47,.92,expression(tau[14]),col='steelblue4')
  }
}
if (sys.nframe() == 0L) {
  dir.create('results/tables', recursive = TRUE, showWarnings = FALSE)
  write.csv(imap_results$summary, 'results/tables/002-imap-cases.csv', row.names = FALSE)
  write.csv(imap_results$errors, 'results/tables/002-imap-errors.csv', row.names = FALSE)
  print(imap_results$summary); print(imap_results$errors)
  writeLines(c(paste('seed:', imap_results$seed), paste('tolerance:', imap_results$tolerance),
               capture.output(sessionInfo())), 'results/002-imap-session.txt')
}
