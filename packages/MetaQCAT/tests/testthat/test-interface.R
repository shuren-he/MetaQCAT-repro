test_that("MetaQCAT exposes a focused bootstrap API", {
  expect_identical(names(formals(MetaQCAT)), c(
    "counts", "X", "X.index", "cluster.id", "taxonomy", "method",
    "min.depth", "min.prev", "n.boot", "fdr.alpha", "save.coef",
    "n.threads", "seed"
  ))
  expect_identical(formals(MetaQCAT)$min.prev, 0.1)
  expect_identical(formals(MetaQCAT)$n.boot, 999L)
  expect_identical(formals(MetaQCAT)$seed, 1234)
  for (obsolete in c("n.perm", "permute.strata", "OTU", "Tax", "cluster",
                     "meta.method", "prev.filter", "save.desc", "nthreads")) {
    a <- list(counts = list(), X = list(), X.index = 1)
    a[[obsolete]] <- TRUE
    expect_error(do.call(MetaQCAT, a), "unused argument",
                 info = obsolete)
  }
})

test_that("advanced interfaces use n.boot and explicit tested indices", {
  d <- make_input()
  x <- lapply(d$X, function(a) cbind(intercept = 1, a))
  for (f in list(Score_test_QCAT_meta)) {
    fit <- f(d$counts, x, 2, n.boot = 9L, save.coef = TRUE)
    expect_true(is.finite(fit$pvalue.comp[["Boot"]]))
  }
  z <- lapply(d$counts, function(y) (y == 0) * 1)
  fit <- Score_test_zero_meta(z, x, 2, n.boot = 9L, save.coef = TRUE)
  expect_true(is.finite(fit$pvalue.pa[["Boot"]]))
  for (f in list(Score_test_QCAT_stat))
    expect_true(is.finite(f(d$counts, x, 2)$score.stat.meta))
  expect_true(is.finite(Score_test_zero_stat(z, x, 2)$score.stat.meta))
  expect_error(Score_test_QCAT_meta(d$counts, x, 1, n.boot = 9L),
               "intercept")
})

test_that("redundant noncontrast count coordinates are not silently generalized", {
  d <- make_input()
  x <- lapply(d$X, function(a) cbind(intercept = 1, a))
  # With a fixed lineage total, all raw taxon-score coordinates obey one
  # sum constraint. The ALR interface removes that redundancy; this legacy
  # noncontrast interface must report it rather than use a pseudoinverse.
  expect_error(Score_test_meta(d$counts, x, 2L, n.boot = 9L),
               "positive definite|singular|rank|estimab")
  expect_error(Score_test_stat(d$counts, x, 2L),
               "positive definite|singular|rank|estimab")
})

test_that("input errors do not impose between/within shuffle assumptions", {
  d <- make_input()
  ids <- lapply(d$X, function(x) factor(rep(1:20, each = 2)))
  fit <- MetaQCAT(d$counts, d$X, 1, cluster.id = ids, n.boot = 9L)
  expect_null(fit$error)
  within.x <- lapply(d$X, function(x) {x[, 1] <- rep(0:1, 20); x})
  fit <- MetaQCAT(d$counts, within.x, 1, cluster.id = ids,
                  n.boot = 9L)
  expect_null(fit$error)

  bad.x <- d$X
  bad.x[[1]] <- bad.x[[1]][40:1, , drop = FALSE]
  expect_error(MetaQCAT(d$counts, bad.x, 1, n.boot = 9L), "row names")
  bad.x <- lapply(d$X, function(x) cbind(1, x))
  expect_error(MetaQCAT(d$counts, bad.x, 2, n.boot = 9L),
               "rank deficient")
  bad.y <- d$counts
  bad.y[[1]][1, 1] <- NA
  expect_error(MetaQCAT(bad.y, d$X, 1, n.boot = 9L), "Invalid")
  expect_error(MetaQCAT(d$counts, d$X, 1, n.boot = 0), "n.boot")
  expect_error(MetaQCAT(d$counts, d$X, 1, n.boot = NA_integer_), "n.boot")
  expect_error(MetaQCAT(d$counts, d$X, 1, min.prev = c(.1, .2),
                        n.boot = 9L), "two values")
})

test_that("one cluster per sample is the default", {
  d <- make_input()
  a <- MetaQCAT(d$counts, d$X, 1, n.boot = 9L, seed = 94L)
  ids <- lapply(d$counts, rownames)
  b <- MetaQCAT(d$counts, d$X, 1, cluster.id = ids,
                n.boot = 9L, seed = 94L)
  expect_equal(a$score.stat.comp, b$score.stat.comp)
  expect_equal(a$bootstrap.stat.comp, b$bootstrap.stat.comp)
  expect_equal(a$pvalue.comp, b$pvalue.comp)
  expect_equal(unname(a$bootstrap$n.clusters), c(40L, 40L))
})
