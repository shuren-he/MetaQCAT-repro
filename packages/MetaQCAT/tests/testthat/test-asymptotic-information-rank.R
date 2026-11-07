test_that("singular cluster covariances are nonestimable for both components", {
  d <- make_input()
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))
  cluster <- lapply(d$counts, function(y) rep("one_subject", nrow(y)))
  z <- lapply(d$counts, function(y) 1L * (y == 0))
  for (method in c("FE-MetaQCAT", "RE-MetaQCAT"))
    expect_error(Score_test_QCAT_stat(d$counts, X, 2L,
                                      cluster.id = cluster, method = method),
                 "positive definite|singular|rank|ill-conditioned|estimab")
  for (method in c("FE-MetaQCAT*", "RE-MetaQCAT*"))
    expect_error(Score_test_zero_stat(z, X, 2L,
                                      cluster.id = cluster, method = method),
                 "positive definite|singular|rank|ill-conditioned|estimab")
})

test_that("unsupported rank does not produce a pseudoinverse p-value", {
  d <- make_input()
  cluster <- lapply(d$counts, function(y) rep("one_subject", nrow(y)))
  expect_warning(
    fit <- MetaQCAT(d$counts, d$X, 1L, taxonomy = d$taxonomy,
                    cluster.id = cluster, method = "FE-MetaQCAT-O", n.boot = 7L),
    "NA columns retained"
  )
  expect_true(all(is.na(fit$lineage.pval)))
  expect_setequal(names(fit$failed.lineage), c("FamilyA", "FamilyB"))
  expect_length(fit$sig.lineage, 0L)
  expect_true(all(vapply(fit$lineage.coef,
                         function(x) x$bootstrap$completed == 0L, logical(1))))
})

test_that("loss of all controls after depth filtering is nonestimable", {
  d <- make_input()
  y <- d$counts
  y[[1L]][1:20, ] <- 0
  expect_warning(
    fit <- MetaQCAT(y, d$X, 1L, n.boot = 7L),
    "Test failed; returning NA"
  )
  expect_identical(fit$status, "nonestimable")
  expect_true(all(is.na(fit$pvalue.comp)))
  expect_equal(fit$bootstrap$completed, 0L)
})
