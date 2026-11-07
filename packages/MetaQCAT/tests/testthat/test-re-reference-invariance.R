test_that("composition RE kernel uses taxon-major coefficient order", {
  q <- 2L
  K <- 3L
  beta <- c(0.3, -0.2, 0.7, 0.1, -0.4, 0.8)
  V <- diag(length(beta))
  parts <- re_qcat_components(beta, V, q)
  expected <- kronecker(diag(K) + matrix(1, K, K), diag(q))
  expect_equal(parts$kernel, expected)
  W <- solve(V)
  z <- solve(V, beta)
  WB <- W %*% expected
  expect_equal(parts$U_theta,
               0.5 * (drop(crossprod(z, expected %*% z)) - sum(diag(WB))))
  expect_equal(parts$V_theta, 0.5 * sum(diag(WB %*% WB)))
})

test_that("RE score transports under a nonorthogonal ALR reference change", {
  set.seed(9303)
  q <- 2L
  beta <- c(0.8, -1.1, 0.5, 0.7, -0.4, 1.3)
  A <- matrix(rnorm(36), 6L, 6L)
  V <- crossprod(A) + diag(0.7, 6L)
  Mtax <- rbind(c(-1, 0, 0), c(-1, 1, 0), c(-1, 0, 1))
  M <- kronecker(Mtax, diag(q))
  beta.new <- drop(M %*% beta)
  V.new <- M %*% V %*% t(M)
  old <- re_qcat_components(beta, V, q)
  new <- re_qcat_components(beta.new, V.new, q)
  expect_equal(new$kernel, M %*% old$kernel %*% t(M), tolerance = 1e-10)
  expect_equal(new$U_theta, old$U_theta, tolerance = 1e-8)
  expect_equal(new$V_theta, old$V_theta, tolerance = 1e-8)
  expect_equal(
    drop(crossprod(beta.new, solve(V.new, beta.new))) +
      re_qcat_statistic(new$U_theta, new$V_theta),
    drop(crossprod(beta, solve(V, beta))) +
      re_qcat_statistic(old$U_theta, old$V_theta),
    tolerance = 1e-8
  )
})

test_that("RE boundary is applied to the summed heterogeneity score", {
  expect_equal(re_qcat_statistic(-2, 3), 0)
  expect_equal(re_qcat_statistic(0, 3), 0)
  expect_equal(re_qcat_statistic(2, 3), 4 / 3)
  expect_equal(re_qcat_statistic(-2 + 1, 3 + 2), 0)
  expect_equal(re_qcat_statistic(-1 + 2, 3 + 2), 1 / 5)
  expect_error(re_qcat_statistic(1, 0),
               "information|positive|degenerate")
})

re_reference_fit <- function(counts, X, tested, reference, n.boot = 19L) {
  .Score_test_QCAT_meta(
    Y.list = counts, X.list = X, X.par.index = tested,
    prev.filter = 0.1, n.boot = n.boot, meta.method = "RE-MetaQCAT",
    seed = 4301L, reference.taxon = reference
  )
}

test_that("all eligible references give the same RE bootstrap p-value", {
  d <- make_input()
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))
  fit <- lapply(colnames(d$counts[[1L]]), function(reference)
    re_reference_fit(d$counts, X, 2L, reference))
  expect_identical(vapply(fit, `[[`, character(1), "reference.taxon"),
                   colnames(d$counts[[1L]]))
  stat <- vapply(fit, `[[`, numeric(1), "score.stat.comp")
  p <- vapply(fit, function(x) unname(x$pvalue.comp[["Boot"]]), numeric(1))
  expect_lt(max(abs(stat - stat[1L])), 1e-6)
  expect_equal(max(abs(p - p[1L])), 0)
})

test_that("every RE bootstrap replicate is reference invariant with q greater than one", {
  d <- make_input()
  X <- lapply(d$X, function(x)
    cbind(intercept = 1, x, trend = seq_len(nrow(x)) %% 3L))
  ids <- lapply(d$counts, function(y) factor(rownames(y)))
  weights <- .make_cluster_weights(ids, n.boot = 19L, seed = 9412L)
  fit <- lapply(c("A", "B", "C"), function(reference) {
    stat <- .Score_test_QCAT_stat(
      Y.list = d$counts, X.list = X, X.par.index = 2:3,
      prev.filter = 0.1, cluster.list = ids,
      meta.method = "RE-MetaQCAT", reference.taxon = reference
    )
    score_bootstrap_meta(
      stat, .weights_for_summary(stat, weights),
      stat$col.index.list, stat$n.par.interest.beta, 2L,
      "RE-MetaQCAT", "count"
    )
  })
  observed <- vapply(fit, `[[`, numeric(1), "observed.statistic")
  expect_lt(max(abs(observed - observed[1L])), 1e-6)
  for (i in 2:3) {
    expect_equal(fit[[i]]$bootstrap.statistic,
                 fit[[1]]$bootstrap.statistic, tolerance = 1e-6)
  }
})

test_that("RE bootstrap stays reference invariant with local missing taxa", {
  d <- make_input()
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))
  counts <- d$counts
  counts[[1]][, "F"] <- 0
  counts[[2]][, "E"] <- 0
  a <- re_reference_fit(counts, X, 2L, "A")
  b <- re_reference_fit(counts, X, 2L, "B")
  expect_equal(a$score.stat.comp, b$score.stat.comp, tolerance = 1e-6)
  expect_equal(a$pvalue.comp, b$pvalue.comp)
})
