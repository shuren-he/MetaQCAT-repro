test_that("one Rademacher draw is shared by every sample of a cluster", {
  ids <- list(
    Study1 = factor(c("a", "a", "b", "b", "c")),
    Study2 = factor(c("a", "b", "b"))
  )
  set.seed(521)
  before <- .Random.seed
  weights <- .make_cluster_weights(ids, n.boot = 7L, seed = 71L)
  expect_identical(.Random.seed, before)
  expect_identical(dim(weights$Study1), c(3L, 7L))
  expect_identical(rownames(weights$Study1), c("a", "b", "c"))
  expect_identical(dim(weights$Study2), c(2L, 7L))
  expect_true(all(unlist(weights) %in% c(-1, 1)))
  sample_weights <- weights$Study1[as.character(ids$Study1), , drop = FALSE]
  expect_identical(sample_weights[1, ], sample_weights[2, ])
  expect_identical(sample_weights[3, ], sample_weights[4, ])
})

test_that("component sample filters keep original study-cluster weight keys", {
  ids <- list(
    Study1 = factor(c("a", "a", "b", "c", "c")),
    Study2 = factor(c("a", "b", "b", "c"))
  )
  weights <- .make_cluster_weights(ids, n.boot = 11L, seed = 311L)
  count_summary <- list(cluster.list = list(
    Study1 = ids$Study1[c(1, 2, 3, 4)],
    Study2 = ids$Study2[c(1, 2, 3)]
  ))
  zero_summary <- list(cluster.list = list(
    Study1 = ids$Study1[c(1, 3, 5)],
    Study2 = ids$Study2[c(2, 4)]
  ))
  cw <- .weights_for_summary(count_summary, weights)
  zw <- .weights_for_summary(zero_summary, weights)
  expect_identical(cw[[1]]["a", ], zw[[1]]["a", ])
  expect_identical(cw[[1]]["b", ], zw[[1]]["b", ])
  expect_identical(cw[[2]]["b", ], zw[[2]]["b", ])
  expect_identical(cw[[1]], weights$Study1)
  expect_identical(rownames(zw[[2]]), c("b", "c"))
})

test_that("fixed-weight study estimates equal hand-multiplied efficient scores", {
  d <- make_input()
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))
  ids <- lapply(d$counts, function(y)
    factor(rep(seq_len(nrow(y) / 2L), each = 2L)))
  weights <- .make_cluster_weights(ids, n.boot = 5L, seed = 109L)
  summary <- .Score_test_QCAT_stat(
    Y.list = d$counts, X.list = X, X.par.index = 2L,
    prev.filter = 0.1, cluster.list = ids,
    meta.method = "FE-MetaQCAT", reference.taxon = "A"
  )
  used_weights <- .weights_for_summary(summary, weights)
  out <- score_bootstrap_meta(
    summary, used_weights, summary$col.index.list,
    summary$n.par.interest.beta, 1L, "FE-MetaQCAT", "count"
  )
  expect_length(out$bootstrap.statistic, 5L)
  expect_equal(out$observed.statistic, summary$score.stat.meta,
               tolerance = 1e-8)
  for (j in seq_along(summary$cluster.score.lst)) {
    U <- summary$cluster.score.lst[[j]]
    bread <- summary$bread.lst[[j]]
    all_weights <- cbind(observed = rep(1, ncol(U)), used_weights[[j]])
    manual <- bread %*% U %*% all_weights
    expect_equal(unname(out$study.beta.reps[[j]]), unname(manual),
                 tolerance = 1e-10)
    expect_equal(drop(manual[, 1L]), drop(summary$est.beta.lst[[j]]),
                 tolerance = 1e-10)
    expect_equal(unname(summary$cov.lst[[j]]),
                 bread %*% tcrossprod(U) %*% t(bread),
                 tolerance = 1e-8)
    expect_equal(summary$psi.lst[[j]], t(bread %*% U), tolerance = 1e-10,
                 ignore_attr = TRUE)
    expect_equal(crossprod(summary$psi.lst[[j]]), summary$cov.lst[[j]],
                 tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("zero-component study estimates use the same cluster score rule", {
  d <- make_input()
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))
  zeros <- lapply(d$counts, function(y) 1L * (y == 0))
  ids <- lapply(d$counts, function(y)
    factor(rep(seq_len(nrow(y) / 2L), each = 2L)))
  weights <- .make_cluster_weights(ids, n.boot = 5L, seed = 81L)
  summary <- .Score_test_zero_stat(
    Y.bin.list = zeros, Z.list = X, Z.par.index = 2L,
    prev.filter = 0.1, cluster.list = ids,
    meta.method = "FE-MetaQCAT*"
  )
  used.weights <- .weights_for_summary(summary, weights)
  out <- score_bootstrap_meta(
    summary, used.weights, summary$col.index.list,
    summary$n.par.interest.alpha, 1L, "FE-MetaQCAT*", "zero"
  )
  expect_length(out$bootstrap.statistic, 5L)
  expect_equal(out$observed.statistic, summary$score.stat.meta,
               tolerance = 1e-8)
  for (j in seq_along(summary$cluster.score.lst)) {
    U <- summary$cluster.score.lst[[j]]
    bread <- summary$bread.lst[[j]]
    manual <- bread %*% U %*%
      cbind(rep(1, ncol(U)), used.weights[[j]])
    expect_equal(unname(out$study.beta.reps[[j]]), unname(manual),
                 tolerance = 1e-10)
    expect_equal(drop(manual[, 1L]), drop(summary$est.alpha.lst[[j]]),
                 tolerance = 1e-10)
    expect_equal(summary$psi.lst[[j]], t(bread %*% U), tolerance = 1e-10,
                 ignore_attr = TRUE)
    expect_equal(crossprod(summary$psi.lst[[j]]), summary$cov.lst[[j]],
                 tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("marginal bootstrap uses every repetition and plus-one correction", {
  d <- make_input()
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))
  for (method in c("FE-MetaQCAT", "RE-MetaQCAT")) {
    fit <- .Score_test_QCAT_meta(
      Y.list = d$counts, X.list = X, X.par.index = 2L,
      n.boot = 13L, meta.method = method, seed = 200L
    )
    expect_length(fit$bootstrap.stat.comp, 13L)
    expected <- (1 + sum(fit$bootstrap.stat.comp >= fit$score.stat.comp)) / 14
    expect_equal(unname(fit$pvalue.comp[["Boot"]]), expected)
    expect_equal(fit$bootstrap$requested, 13L)
    expect_equal(fit$bootstrap$completed, 13L)
  }
})

test_that("ties are included and an all-one plan runs every requested draw", {
  d <- make_input()
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))
  ids <- .complete_cluster_ids(d$counts, NULL)
  weights <- lapply(ids, function(id)
    matrix(1, nlevels(id), 7L, dimnames = list(levels(id), NULL)))
  fit <- .Score_test_QCAT_meta(
    Y.list = d$counts, X.list = X, X.par.index = 2L,
    cluster.list = ids, weights = weights, n.boot = 7L,
    meta.method = "FE-MetaQCAT", seed = 71L
  )
  expect_length(fit$bootstrap.stat.comp, 7L)
  expect_equal(fit$bootstrap.stat.comp, rep(fit$score.stat.comp, 7L),
               tolerance = 1e-10)
  expect_equal(fit$pvalue.comp[["Boot"]], 1)
  expect_equal(fit$bootstrap$completed, 7L)
})

test_that("omnibus combines marginal plus-one p-values and shares weighted draws", {
  d <- make_input()
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))
  ids <- lapply(d$counts, function(y)
    factor(rep(seq_len(nrow(y) / 2L), each = 2L)))
  weights <- .make_cluster_weights(ids, n.boot = 7L, seed = 317L)
  zeros <- lapply(d$counts, function(y) 1L * (y == 0))
  combined <- .run_meta_test(
    counts = d$counts, zeros = zeros, X = X, X.index = 2L,
    cluster.id = ids, weights = weights, method = "FE-MetaQCAT-O",
    min.prev = c(.1, .1), n.boot = 7L,
    save.coef = FALSE, seed = 317L
  )
  positive <- .Score_test_QCAT_meta(
    Y.list = d$counts, X.list = X, X.par.index = 2L,
    cluster.list = ids, weights = weights, n.boot = 7L,
    meta.method = "FE-MetaQCAT", seed = 317L
  )
  zero <- .Score_test_zero_meta(
    Y.bin.list = zeros, Z.list = X, Z.par.index = 2L,
    cluster.list = ids, weights = weights, n.boot = 7L,
    meta.method = "FE-MetaQCAT*", seed = 317L
  )
  expect_null(combined$error)
  expect_equal(combined$bootstrap.stat.comp, positive$bootstrap.stat.comp)
  expect_equal(combined$bootstrap.stat.pa, zero$bootstrap.stat.pa)
  expected <- ACAT(c(positive$pvalue.comp[["Boot"]],
                     zero$pvalue.pa[["Boot"]]))
  expect_equal(combined$pvalue.comb[["Boot"]], expected)
  expect_null(combined$bootstrap.stat.comb)
  expect_null(combined$bootstrap.omnibus$internal.p)
})

test_that("all sign patterns have exact mean zero and stored covariance", {
  Psi <- cbind(c(.2, -.3, .4, .7, -.5, .1),
                c(.5, .1, -.2, .3, .6, -.4))
  signs <- t(as.matrix(expand.grid(rep(list(c(-1, 1)), nrow(Psi)))))
  summary <- list(cluster.score.lst = list(t(Psi)),
                  bread.lst = list(diag(ncol(Psi))),
                  cov.lst = list(crossprod(Psi)))
  out <- score_bootstrap_meta(summary, list(signs), list(1:2),
                              2L, 1L, "FE-MetaQCAT", "count")
  beta <- out$study.beta.reps[[1L]][, -1L, drop = FALSE]
  expect_equal(unname(beta), unname(crossprod(Psi, signs)), tolerance = 1e-10)
  expect_equal(rowMeans(beta), c(0, 0), tolerance = 1e-10)
  # Population normalization is required when enumerating the full sign law.
  expect_equal(tcrossprod(beta) / ncol(beta), crossprod(Psi), tolerance = 1e-10)
  expect_equal(drop(out$study.beta.reps[[1L]][, 1L]), colSums(Psi),
               tolerance = 1e-10)
  all.one <- which(colSums(signs) == nrow(signs))
  expect_length(all.one, 1L)
  expect_equal(out$bootstrap.statistic[all.one], out$observed.statistic,
               tolerance = 1e-10)
})

test_that("nonfinite bootstrap draws fail instead of changing the denominator", {
  Psi <- cbind(c(.2, -.3, .4, .7, -.5, .1),
                c(.5, .1, -.2, .3, .6, -.4))
  summary <- list(cluster.score.lst = list(t(Psi)),
                  bread.lst = list(diag(ncol(Psi))),
                  cov.lst = list(crossprod(Psi)))
  weights <- matrix(1, nrow(Psi), 3L)
  weights[2L, 2L] <- NA_real_
  expect_error(score_bootstrap_meta(summary, list(weights), list(1:2),
                                     2L, 1L, "FE-MetaQCAT", "count"),
               "non-finite|Non-finite")
  weights[2L, 2L] <- 0
  expect_error(score_bootstrap_meta(summary, list(weights), list(1:2),
                                     2L, 1L, "FE-MetaQCAT", "count"),
               "Rademacher")
})

test_that("count and zero RE use the positive part for every observed and bootstrap score", {
  for (component in c("count", "zero")) {
    suffix <- if (component == "zero") "*" else ""
    seen <- c(negative = FALSE, zero = FALSE, positive = FALSE)
    for (psi in list(c(1, -1, .5, -.5, .25, -.25),
                     c(1, 0, 0, 0, 0, 0), rep(.5, 6L))) {
      Psi <- matrix(psi, ncol = 1L)
      summary <- list(cluster.score.lst = list(t(Psi)),
                      bread.lst = list(matrix(1)),
                      cov.lst = list(crossprod(Psi)))
      weights <- t(as.matrix(expand.grid(rep(list(c(-1, 1)), 6L))))
      fe <- score_bootstrap_meta(summary, list(weights), list(1L),
                                 1L, 1L, paste0("FE-MetaQCAT", suffix), component)
      re <- score_bootstrap_meta(summary, list(weights), list(1L),
                                 1L, 1L, paste0("RE-MetaQCAT", suffix), component)
      u <- drop(re$re.score)
      seen <- seen | c(any(u < -1e-12), any(abs(u) <= 1e-12), any(u > 1e-12))
      increment <- pmax(0, u)^2 / re$re.information
      expect_equal(c(re$observed.statistic, re$bootstrap.statistic),
                   c(fe$observed.statistic, fe$bootstrap.statistic) + increment,
                   tolerance = 1e-10, info = component)
      expect_length(re$bootstrap.statistic, ncol(weights))
    }
    expect_true(all(seen), info = component)
  }
})

test_that("both RE components truncate after summing opposing study scores", {
  psi <- list(rep(.5, 6L), c(.1, -.1, .05, -.05, .025, -.025))
  summary <- list(cluster.score.lst = lapply(psi, function(x) matrix(x, nrow = 1L)),
                  bread.lst = list(matrix(1), matrix(1)),
                  cov.lst = lapply(psi, function(x) matrix(sum(x^2))))
  weights <- t(as.matrix(expand.grid(rep(list(c(-1, 1)), 6L))))
  beta <- lapply(psi, function(x) c(sum(x), drop(x %*% weights)))
  variance <- vapply(psi, function(x) sum(x^2), numeric(1))
  for (component in c("count", "zero")) {
    suffix <- if (component == "zero") "*" else ""
    kernel <- if (component == "count") 2 else 1
    u <- lapply(seq_along(psi), function(i)
      .5 * kernel * (beta[[i]]^2 / variance[i]^2 - 1 / variance[i]))
    total.u <- Reduce(`+`, u)
    total.v <- .5 * kernel^2 * sum(1 / variance^2)
    expect_lt(total.u[1L], 0)
    expect_gt(sum(vapply(u, function(x) max(0, x[1L]), numeric(1))), 0)
    fe <- score_bootstrap_meta(summary, list(weights, weights), list(1L, 1L),
                               1L, 1L, paste0("FE-MetaQCAT", suffix), component)
    re <- score_bootstrap_meta(summary, list(weights, weights), list(1L, 1L),
                               1L, 1L, paste0("RE-MetaQCAT", suffix), component)
    expect_equal(drop(re$re.score), total.u, tolerance = 1e-10)
    expect_equal(re$re.information, total.v, tolerance = 1e-10)
    expect_equal(c(re$observed.statistic, re$bootstrap.statistic),
                 c(fe$observed.statistic, fe$bootstrap.statistic) +
                   pmax(0, total.u)^2 / total.v, tolerance = 1e-10)
    expect_equal(re$observed.statistic, fe$observed.statistic)
  }
})

test_that("zero RE diagnostics record the one-sided boundary convention", {
  d <- make_input()
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))
  z <- lapply(d$counts, function(y) 1L * (y == 0))
  out <- .Score_test_zero_stat(z, X, 2L, meta.method = "RE-MetaQCAT*")
  expect_identical(out$re.diagnostics$kernel, "identity")
  expect_identical(out$re.diagnostics$variance.score.formula, "positive-part")
  expect_identical(out$re.diagnostics$positive.part.active,
                   out$re.diagnostics$U_theta > 0)
})

test_that("the separate legacy Poisson helper is unchanged", {
  Psi <- matrix(c(1, -1, .5, -.5, .25, -.25), ncol = 1L)
  summary <- list(cluster.score.lst = list(t(Psi)), bread.lst = list(matrix(1)),
                  cov.lst = list(crossprod(Psi)))
  weights <- t(as.matrix(expand.grid(rep(list(c(-1, 1)), 6L))))
  fe <- score_bootstrap_meta(summary, list(weights), list(1L),
                             1L, 1L, "FE-MetaQCAT", "poisson")
  re <- score_bootstrap_meta(summary, list(weights), list(1L),
                             1L, 1L, "RE-MetaQCAT", "poisson")
  expect_lt(drop(re$re.score)[1L], 0)
  expect_equal(c(re$observed.statistic, re$bootstrap.statistic),
               c(fe$observed.statistic, fe$bootstrap.statistic) +
                 drop(re$re.score)^2 / re$re.information, tolerance = 1e-10)
})
