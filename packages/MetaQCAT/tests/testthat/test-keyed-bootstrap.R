test_that("study and subject order do not reassign keyed multipliers", {
  ids <- list(
    StudyZ = factor(c("c", "a", "a", "b"), levels = c("c", "a", "b")),
    StudyA = factor(c("y", "x", "x"), levels = c("y", "x"))
  )
  a <- .make_cluster_weights(ids, n.boot = 13L, seed = 731L)
  reordered <- rev(lapply(ids, function(id)
    factor(rev(as.character(id)), levels = rev(levels(id)))))
  b <- .make_cluster_weights(reordered, n.boot = 13L, seed = 731L)
  for (study in names(a))
    expect_equal(a[[study]], b[[study]][rownames(a[[study]]), , drop = FALSE])
})

test_that("independent observations use stable sample keys", {
  d <- make_input()
  a <- .complete_cluster_ids(d$counts, NULL)
  b <- .complete_cluster_ids(rev(lapply(d$counts, function(y)
    y[nrow(y):1L, , drop = FALSE])), NULL)
  for (study in names(a)) {
    expect_setequal(as.character(a[[study]]), rownames(d$counts[[study]]))
    expect_equal(rev(as.character(b[[study]])), as.character(a[[study]]))
  }
  wa <- .make_cluster_weights(a, n.boot = 11L, seed = 735L)
  wb <- .make_cluster_weights(b, n.boot = 11L, seed = 735L)
  for (study in names(wa))
    expect_equal(wa[[study]], wb[[study]][rownames(wa[[study]]), , drop = FALSE])
})

test_that("fixed preprocessed outcomes give the same fit after visit reordering", {
  d <- make_input()
  x <- lapply(d$X, function(a) cbind(intercept = 1, a))
  ids <- lapply(d$counts, function(y)
    factor(paste0("subject", rep(seq_len(nrow(y) / 2L), each = 2L))))
  z <- lapply(d$counts, function(y) 1L * (y == 0))
  order <- c(seq(2L, 40L, by = 2L), seq(1L, 39L, by = 2L))
  reorder <- function(a) rev(lapply(a, function(b) {
    if (is.matrix(b)) b[order, , drop = FALSE] else factor(b[order])
  }))
  run <- function(y, zero, design, id, method) {
    weights <- .make_cluster_weights(id, n.boot = 19L, seed = 912L)
    .run_meta_test(counts = y, zeros = zero, X = design, X.index = 2L,
                   cluster.id = id, weights = weights, method = method,
                   min.prev = c(.1, .1), n.boot = 19L,
                   save.coef = FALSE, seed = 912L)
  }
  for (method in c("FE-MetaQCAT-O", "RE-MetaQCAT-O")) {
    a <- run(d$counts, z, x, ids, method)
    b <- run(reorder(d$counts), reorder(z), reorder(x), reorder(ids), method)
    expect_null(a$error)
    expect_null(b$error)
    for (field in c("score.stat.comp", "score.stat.pa",
                     "bootstrap.stat.comp", "bootstrap.stat.pa"))
      expect_equal(a[[field]], b[[field]], tolerance = 1e-7, info = field)
    expect_equal(a$pvalue.comp, b$pvalue.comp)
    expect_equal(a$pvalue.pa, b$pvalue.pa)
    expect_equal(a$pvalue.comb, b$pvalue.comb)
  }
})

test_that("missing and empty subject IDs fail explicitly", {
  d <- make_input()
  id <- lapply(d$counts, rownames)
  id[[1L]][2L] <- NA_character_
  expect_error(MetaQCAT(d$counts, d$X, 1L, cluster.id = id, n.boot = 5L),
               "ID|id")
  id[[1L]][2L] <- ""
  expect_error(MetaQCAT(d$counts, d$X, 1L, cluster.id = id, n.boot = 5L),
               "ID|id")
  y <- d$counts
  x <- d$X
  rownames(y[[1L]])[2L] <- rownames(y[[1L]])[1L]
  rownames(x[[1L]]) <- rownames(y[[1L]])
  expect_error(MetaQCAT(y, x, 1L, n.boot = 5L),
               "[Ss]ample|[Uu]nique|duplicat")
})

test_that("canonical preprocessing preserves keyed results after all input reordering", {
  d <- make_input()
  ids <- lapply(d$counts, function(y)
    paste0("subject", rep(seq_len(nrow(y) / 2L), each = 2L)))
  fit <- MetaQCAT(d$counts, d$X, 1L, taxonomy = d$taxonomy,
                  cluster.id = ids, method = "RE-MetaQCAT-O",
                  n.boot = 13L, seed = 327L)
  y <- rev(lapply(d$counts, function(a) a[40:1, 6:1, drop = FALSE]))
  x <- rev(lapply(d$X, function(a) a[40:1, , drop = FALSE]))
  ids <- rev(lapply(ids, rev))
  reordered <- MetaQCAT(y, x, 1L, taxonomy = d$taxonomy[6:1, , drop = FALSE],
                        cluster.id = ids, method = "RE-MetaQCAT-O",
                        n.boot = 13L, seed = 327L)
  expect_length(fit$failed.lineage, 0L)
  expect_length(reordered$failed.lineage, 0L)
  expect_equal(fit$lineage.pval, reordered$lineage.pval)
  for (family in names(fit$lineage.coef)) {
    a <- fit$lineage.coef[[family]]
    b <- reordered$lineage.coef[[family]]
    expect_equal(a$bootstrap.stat.comp, b$bootstrap.stat.comp, tolerance = 1e-8)
    expect_equal(a$bootstrap.stat.pa, b$bootstrap.stat.pa, tolerance = 1e-8)
  }
})
