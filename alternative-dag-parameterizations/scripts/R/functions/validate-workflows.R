imap_validate_frugal <- function(seed = 20261011L, cases = 100L, tolerance = 1e-9) {
  set.seed(seed)
  # The last three cases retain legitimate zero-edge parameter values.
  coefficients <- rbind(matrix(runif(2 * cases, -1, 1), cases, 2),
                        c(0, .8), c(.8, 0), c(0, 0))
  errors <- c(frugal_covariance = 0, inverse = 0, chain_constraint = 0,
              intervention_variance = 0, relabeled_imap = 0, valid_adjustment = 0)
  for (k in seq_len(nrow(coefficients))) {
    b <- coefficients[k, 1]; c <- coefficients[k, 2]; v <- runif(3, .5, 2)
    L <- diag(3); L[1,2] <- -b; L[2,3] <- -c
    z <- imap_chain_frugal(b, c, v); p <- imap_frugal_parameters(z)
    inverse <- imap_frugal_chain_inverse(z)
    errors['frugal_covariance'] <- max(errors['frugal_covariance'],
                                       abs(imap_covariance(L,v) - imap_covariance(p$L,p$v)))
    errors['inverse'] <- max(errors['inverse'], abs(inverse$b-b), abs(inverse$c-c), abs(inverse$v-v))
    errors['chain_constraint'] <- max(errors['chain_constraint'], abs(p$L[1,3]))
    causal <- imap_intervention(L,v,1,3)
    errors['intervention_variance'] <- max(errors['intervention_variance'],
                                           abs(causal$effect-z$tau), abs(causal$variance-z$s))
    altered <- z; altered$tau <- z$tau + .2
    rejected <- tryCatch(imap_frugal_chain_inverse(altered), error = function(e) e)
    stopifnot(inherits(rejected, 'error'))
  }
  # Every relabeling of the diamond checks arbitrary input topological orders.
  permutations <- function(x) {
    if (length(x) == 1L) return(matrix(x, 1L))
    do.call(rbind, lapply(x, function(j) cbind(j, permutations(setdiff(x, j)))))
  }
  perm <- permutations(1:4); A <- matrix(FALSE,4,4)
  A[rbind(c(1,2), c(1,3), c(2,4), c(3,4))] <- TRUE
  L <- diag(4); L[A] <- c(-.7,.4,-.6,.8); v <- c(1.2,.7,1.5,.9)
  for (k in seq_len(nrow(perm))) {
    label <- perm[k,]; An <- A; Ln <- L; vn <- v
    An[label,label] <- A; Ln[label,label] <- L; vn[label] <- v
    p <- imap_transform(Ln,vn,An,label[1],label[4])
    reference <- imap_cholesky(imap_covariance(Ln,vn), p$ordering)
    errors['relabeled_imap'] <- max(errors['relabeled_imap'], abs(p$L-reference$L),
      abs(p$v-reference$v), abs(-p$L[label[1],label[4]]-solve(L)[1,4]))
    o <- imap_orders(An,label[1],label[4]); G <- imap_ordered_graph(An,o$optimal)
    stopifnot(setequal(which(G[,label[4]]), c(label[1],o$O)))
  }
  # Candidate-set comparison adapted from the archived two-step example.
  candidates <- list(integer(), 1L, 2L, c(1L,2L), c(1L,3L), c(1L,2L,3L), 5L)
  graphs <- imap_adjustment_example_graphs()
  validity <- matrix(FALSE, length(candidates), length(graphs))
  valid_checks <- 0L
  for (k in seq_along(graphs)) {
    A <- graphs[[k]]; L <- diag(6) - .7 * A; Sigma <- imap_covariance(L,rep(1,6))
    effect <- imap_intervention(L,rep(1,6),4,6)$effect
    for (j in seq_along(candidates)) {
      Z <- candidates[[j]]; validity[j,k] <- imap_valid_pretreatment_adjustment(A,4,6,Z)
      if (validity[j,k]) {
        nodes <- c(4,Z); slope <- solve(Sigma[nodes,nodes,drop=FALSE],Sigma[nodes,6])[1]
        errors['valid_adjustment'] <- max(errors['valid_adjustment'], abs(slope-effect))
        valid_checks <- valid_checks + 1L
      }
    }
  }
  stopifnot(all(errors < tolerance), all(validity[2,]), !any(validity[7,]))
  list(summary = data.frame(check = c('Random chain values', 'Zero-edge chain values',
         'Diamond relabelings', 'Graph/set pairs', 'Valid population regressions'),
         cases = c(cases,3,nrow(perm),length(validity),valid_checks)),
       errors = data.frame(check = names(errors), maximum_absolute_error = unname(errors)),
       validity = data.frame(set = c('empty','1','2','1,2','1,3','1,2,3','5'),
                             setNames(as.data.frame(validity),names(graphs))),
       seed = seed, tolerance = tolerance)
}

imap_validate_bayesian <- function(seed = 20261012L, S = 8000L, n = 120L, tolerance = 1e-9) {
  set.seed(seed)
  # A deliberately small equivalence class; this is enumeration, not graph MCMC.
  chain <- fork <- reverse <- matrix(FALSE,3,3)
  chain[1,2] <- chain[2,3] <- TRUE
  fork[2,1] <- fork[2,3] <- TRUE
  reverse[3,2] <- reverse[2,1] <- TRUE
  graphs <- list(chain = chain, fork = fork, reverse = reverse)
  L <- diag(3); L[1,2] <- -.7; L[2,3] <- -.9
  X <- matrix(rnorm(n*3),n,3) %*% chol(imap_covariance(L,rep(1,3)))
  a <- 8; U <- diag(3); alpha <- dagwishart_shapes(chain,a)
  scores <- vapply(graphs,function(A) dagwishart_log_marginal(X,A,a,U),0.0)
  equivalent_score_error <- max(abs(scores-scores[1]))
  prior <- dagwishart_blocks(chain,alpha,U)
  posterior <- dagwishart_blocks(chain,alpha+n,U+crossprod(X))
  # Separate node blocks are independent, so moments of B12*B23 factor.
  prior_moments <- dagwishart_coefficient_moments(prior[[2]],1) *
                   dagwishart_coefficient_moments(prior[[3]],2)
  posterior_moments <- dagwishart_coefficient_moments(posterior[[2]],1) *
                       dagwishart_coefficient_moments(posterior[[3]],2)
  C <- unname(prior_moments['second'])
  correction <- unname(posterior_moments['second']) / C
  new_scores <- scores; new_scores[1] <- scores[1] + log(correction)
  baseline_mass <- imap_normalized_weights(scores)
  tilted_mass <- imap_normalized_weights(new_scores)
  draws <- dagwishart_sample(S,chain,alpha+n,U+crossprod(X))
  tau <- draws$L[1,2,] * draws$L[2,3,]
  # Check the transformed effect against structural surgery for sample draws.
  transformation_error <- 0
  for (i in seq_len(min(S,100L))) {
    p <- imap_transform(draws$L[,,i],draws$v[i,],chain,1,3)
    transformation_error <- max(transformation_error,abs(-p$L[1,3]-tau[i]),
      abs(imap_intervention(draws$L[,,i],draws$v[i,],1,3)$effect-tau[i]))
  }
  weights <- imap_normalized_weights(2*log(abs(tau)))
  tilted_mean <- unname(posterior_moments['third']/posterior_moments['second'])
  weighted_mean <- sum(weights*tau)
  second_mc <- mean(tau^2); second_se <- sd(tau^2)/sqrt(S)
  # Delta-method Monte Carlo error for the self-normalized tilted mean.
  mean_se <- sd(tau^2*(tau-tilted_mean))/sqrt(S)/mean(tau^2)
  discrepancies <- c(abs(second_mc-posterior_moments['second']), abs(weighted_mean-tilted_mean))
  ses <- c(second_se,mean_se)
  stopifnot(equivalent_score_error < tolerance, transformation_error < tolerance,
            all(discrepancies < 7*ses))
  list(scores = data.frame(graph = names(graphs), baseline_log_evidence = unname(scores),
         baseline_probability = baseline_mass, tilted_probability = tilted_mass),
       moments = data.frame(quantity = c('Prior effect second moment', 'Posterior baseline second moment',
         'Tilted/baseline evidence ratio', 'Active tilted posterior mean', 'Averaged tilted posterior mean'),
         value = c(C,posterior_moments['second'],correction,tilted_mean,tilted_mass[1]*tilted_mean)),
       monte_carlo = data.frame(quantity = c('Posterior effect second moment','Active tilted posterior mean'),
         exact = c(posterior_moments['second'],tilted_mean), estimate = c(second_mc,weighted_mean),
         monte_carlo_se = ses),
       diagnostics = data.frame(check = c('Equivalent-score absolute discrepancy',
         'Effect-transformation absolute discrepancy','Importance effective sample size'),
         value = c(equivalent_score_error,transformation_error,1/sum(weights^2))),
       seed = seed, samples = S, n = n, tolerance = tolerance)
}
