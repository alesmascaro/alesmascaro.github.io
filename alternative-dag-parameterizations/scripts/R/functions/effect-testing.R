# Exact enumeration and DAG-Wishart effect moments for small, fully observed DAGs.
imap_enumerate_dags <- function(q = 3L) {
  stopifnot(q >= 2L, q <= 4L)
  pairs <- t(combn(seq_len(q), 2)); graphs <- list()
  for (code in 0:(3^nrow(pairs) - 1L)) {
    A <- matrix(FALSE, q, q); digits <- code
    for (k in seq_len(nrow(pairs))) {
      state <- digits %% 3L; digits <- digits %/% 3L
      if (state == 1L) A[pairs[k, 1], pairs[k, 2]] <- TRUE
      if (state == 2L) A[pairs[k, 2], pairs[k, 1]] <- TRUE
    }
    R <- imap_closure(A)
    if (!any(diag(R))) graphs[[length(graphs) + 1L]] <- A
  }
  graphs
}

imap_directed_paths <- function(A, d, y) {
  imap_topological_order(A)
  walk <- function(u, nodes) {
    if (u == y) return(list(nodes))
    out <- list()
    for (v in which(A[u, ])) out <- c(out, walk(v, c(nodes, v)))
    out
  }
  walk(d, d)
}

dagwishart_effect_second_moment <- function(blocks, paths) {
  # Node blocks are independent; a pair of paths uses at most two entries per node.
  if (!length(paths)) return(0)
  total <- 0
  for (p in paths) for (r in paths) {
    edges <- rbind(cbind(head(p, -1L), tail(p, -1L)),
                   cbind(head(r, -1L), tail(r, -1L)))
    term <- 1
    for (j in unique(edges[, 2])) {
      z <- blocks[[j]]; k <- match(edges[edges[, 2] == j, 1], z$P)
      mu <- -z$mu[k]
      if (length(k) == 1L) term <- term * mu
      else {
        if (z$r <= 1) stop("Effect second moment needs finite node variances")
        term <- term * (prod(mu) + z$rate / (z$r - 1) * z$Kinv[k[1], k[2]])
      }
    }
    total <- total + term
  }
  if (total < -1e-10) stop("Negative effect second moment")
  max(0, total)
}

imap_effect_testing_setup <- function(a = 8, U = diag(3)) {
  graphs <- imap_enumerate_dags(3L)
  paths <- lapply(graphs, imap_directed_paths, d = 1L, y = 3L)
  active <- lengths(paths) > 0L
  blocks <- lapply(graphs, function(A) dagwishart_blocks(A, dagwishart_shapes(A, a), U))
  constants <- vapply(seq_along(graphs), function(k)
    dagwishart_effect_second_moment(blocks[[k]], paths[[k]]), 0.0)
  stopifnot(length(graphs) == 25L, all(constants[active] > 0), all(constants[!active] == 0))
  list(graphs = graphs, paths = paths, active = active, constants = constants,
       log_prior_normalizers = vapply(graphs, function(A)
         dagwishart_log_normalizer(A, dagwishart_shapes(A, a), U), 0.0), a = a, U = U)
}

imap_effect_testing_scores <- function(X, setup) {
  n <- nrow(X); Upost <- setup$U + crossprod(X)
  log_m <- second <- numeric(length(setup$graphs))
  for (k in seq_along(setup$graphs)) {
    A <- setup$graphs[[k]]; alpha <- dagwishart_shapes(A, setup$a)
    log_m[k] <- -n * ncol(X) / 2 * log(2*pi) +
      dagwishart_log_normalizer(A, alpha + n, Upost) - setup$log_prior_normalizers[k]
    if (setup$active[k]) second[k] <- dagwishart_effect_second_moment(
      dagwishart_blocks(A, alpha + n, Upost), setup$paths[[k]])
  }
  tilted <- log_m
  tilted[setup$active] <- log_m[setup$active] +
    log(second[setup$active] / setup$constants[setup$active])
  # Equal graph mass for both priors: the only difference is the active parameter prior.
  p0 <- imap_normalized_weights(log_m); pM <- imap_normalized_weights(tilted)
  list(baseline = sum(p0[setup$active]), nonlocal = sum(pM[setup$active]),
       log_evidence = log_m, tilted_log_evidence = tilted, second = second)
}

imap_validate_effect_moments <- function(seed = 20261014L, samples = 30000L) {
  set.seed(seed); setup <- imap_effect_testing_setup()
  L <- diag(3); L[1,2] <- L[2,3] <- -.8
  X <- matrix(rnorm(150),50,3) %*% chol(imap_covariance(L,rep(1,3)))
  checks <- list()
  for (k in which(setup$active)) for (posterior in c(FALSE,TRUE)) {
    A <- setup$graphs[[k]]
    alpha <- dagwishart_shapes(A,setup$a) + if(posterior) nrow(X) else 0
    U <- setup$U + if(posterior) crossprod(X) else 0
    draws <- dagwishart_sample(samples,A,alpha,U)
    # Independent three-node path formula, including the direct/indirect cross term.
    tau <- -draws$L[1,3,] + draws$L[1,2,]*draws$L[2,3,]
    exact <- dagwishart_effect_second_moment(dagwishart_blocks(A,alpha,U),setup$paths[[k]])
    se <- sd(tau^2)/sqrt(samples); standardized <- abs(mean(tau^2)-exact)/se
    stopifnot(standardized < 7)
    checks[[length(checks)+1L]] <- data.frame(graph=k,posterior=posterior,
      exact=exact,estimate=mean(tau^2),mc_se=se,standardized_error=standardized)
  }
  list(checks=do.call(rbind,checks), seed=seed, samples=samples)
}

imap_simulate_effect_testing <- function(seed = 20261013L, replicates = 200L,
    sizes = c(25L, 50L, 100L, 250L, 500L, 1000L), a = 8) {
  stopifnot(replicates >= 2L, all(diff(sizes) > 0), min(sizes) > 3)
  set.seed(seed); setup <- imap_effect_testing_setup(a)
  Bchain <- Bcollider <- matrix(0, 3, 3)
  Bchain[1,2] <- Bchain[2,3] <- .8
  Bcollider[1,2] <- Bcollider[3,2] <- .8
  truths <- list(chain = Bchain, collider = Bcollider)
  rows <- vector("list", length(truths) * replicates * length(sizes)); pos <- 0L
  for (scenario in names(truths)) {
    Sigma <- imap_covariance(diag(3) - truths[[scenario]], rep(1,3))
    for (r in seq_len(replicates)) {
      # Nested sizes within a replicate; independent replicates and scenarios.
      X <- matrix(rnorm(max(sizes)*3), max(sizes), 3) %*% chol(Sigma)
      for (n in sizes) {
        z <- imap_effect_testing_scores(X[seq_len(n),,drop=FALSE], setup)
        pos <- pos + 1L
        rows[[pos]] <- data.frame(scenario = scenario, replicate = r, n = n,
                                 baseline = z$baseline, nonlocal = z$nonlocal)
      }
    }
  }
  raw <- do.call(rbind, rows); summary <- list(); pos <- 0L
  for (scenario in names(truths)) for (n in sizes) for (prior in c("baseline", "nonlocal")) {
    p <- raw[raw$scenario == scenario & raw$n == n, prior]
    pos <- pos + 1L
    summary[[pos]] <- data.frame(scenario = scenario, n = n, prior = prior,
      mean = mean(p), mc_se = sd(p)/sqrt(replicates),
      q10 = unname(quantile(p, .1)), q90 = unname(quantile(p, .9)))
  }
  stopifnot(all(is.finite(as.matrix(raw[,c("baseline","nonlocal")]))),
            all(raw$baseline >= 0 & raw$baseline <= 1),
            all(raw$nonlocal >= 0 & raw$nonlocal <= 1))
  list(raw = raw, summary = do.call(rbind, summary), seed = seed,
       replicates = replicates, sizes = sizes, a = a, U = setup$U,
       graph_count = length(setup$graphs), active_count = sum(setup$active),
       prior_active_probability = mean(setup$active),
       truths = truths, error_variances = rep(1,3))
}

imap_plot_effect_testing <- function(result) {
  old <- par(mfrow=c(2,2), mar=c(3.2,3.8,2.5,.7), oma=c(2,0,0,0))
  on.exit(par(old)); colors <- c(baseline="#27649B", nonlocal="#B54831")
  for (event in c("active", "null")) for (scenario in c("chain", "collider")) {
    plot(range(result$sizes), c(0,1), type="n", log="x", xaxt="n",
      xlab="Sample size", ylab=if(event == "active") "Posterior P(effect present)" else "Posterior P(effect absent)",
      main=if(scenario == "chain") "True chain: effect = 0.64" else "True collider: effect = 0")
    axis(1, at=result$sizes, labels=result$sizes, cex.axis=.85)
    abline(h=if(event == "active") result$prior_active_probability else 1-result$prior_active_probability,
           lty=3, col="grey50")
    for (prior in names(colors)) {
      z <- result$summary[result$summary$scenario == scenario & result$summary$prior == prior,]
      center <- z$mean; lo <- z$q10; hi <- z$q90
      if(event == "null") { center <- 1-center; lo <- 1-z$q90; hi <- 1-z$q10 }
      polygon(c(z$n,rev(z$n)), c(lo,rev(hi)), border=NA,
              col=adjustcolor(colors[prior], alpha.f=.13))
      lines(z$n,center,col=colors[prior],lwd=2.2,lty=if(prior=="baseline") 1 else 2)
      points(z$n,center,col=colors[prior],pch=if(prior=="baseline") 16 else 17,cex=.65)
    }
    if(event == "active" && scenario == "chain") legend("bottomleft",
      legend=c("DAG-Wishart", "Normalized non-local tilt"), col=colors, lty=c(1,2), lwd=2,
      bty="n", cex=.8)
  }
  mtext("Lines: replicate means; shading: 10th-90th percentiles; dotted: common prior probability",
        side=1, outer=TRUE, line=.5, cex=.78)
}
