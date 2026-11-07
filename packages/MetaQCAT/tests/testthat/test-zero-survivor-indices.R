test_that("zero-part coordinates come only from successfully fitted studies", {
  d <- make_input()
  zeros <- lapply(d$counts, function(y) (y == 0) * 1)
  zeros[[1L]][, "B"] <- 1
  X <- lapply(d$X, function(x) cbind(intercept = 1, x))

  fit.calls <- 0L
  local_mocked_bindings(
    glm = function(...) {
      fit.calls <<- fit.calls + 1L
      if (fit.calls == 2L) return(list(converged = FALSE))
      stats::glm(...)
    },
    .package = "MetaQCAT"
  )

  expect_warning(
    fit <- Score_test_zero_stat(zeros, X, 2L, min.prev = 0),
    "Cannot converge for study 2"
  )
  expect_equal(fit.calls, 2L)
  expect_length(fit$Y.bin.list, 1L)
  expect_identical(colnames(fit$alpha.hat), setdiff(colnames(zeros[[1L]]), "B"))
  expect_identical(fit$n.par.interest.alpha, 5L)
  expect_identical(as.integer(fit$col.index.list[[1L]]), seq_len(5L))
})
