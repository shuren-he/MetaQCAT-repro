test_that("six methods return the requested fixed-budget bootstrap", {
  d <- make_input()
  methods <- c("FE-MetaQCAT", "FE-MetaQCAT*", "FE-MetaQCAT-O",
               "RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O")
  for (m in methods) {
    fit <- MetaQCAT(d$counts, d$X, 1, method = m, taxonomy = d$taxonomy,
                    n.boot = 9L, save.coef = TRUE, seed = 71L)
    expect_length(fit$failed.lineage, 0L)
    expect_identical(colnames(fit$lineage.pval), c("FamilyA", "FamilyB"))
    expect_true(all(is.finite(fit$lineage.pval)), info = m)
    expect_true(all(fit$lineage.pval >= 0 & fit$lineage.pval <= 1), info = m)
    expect_true(any(grepl("^Bootstrap-", rownames(fit$lineage.pval))),
                info = m)
    expect_identical(dim(fit$lineage.qval), dim(fit$lineage.pval))
    expect_equal(fit$bootstrap$requested, 9L, info = m)
    expect_equal(fit$bootstrap$seed, 71L, info = m)

    global <- MetaQCAT(d$counts, d$X, 1, method = m, n.boot = 9L,
                       seed = 71L)
    expect_null(global$error, info = m)
    p <- if (grepl("-O$", m)) global$pvalue.comb else if (grepl("\\*$", m))
      global$pvalue.pa else global$pvalue.comp
    expect_true(is.finite(p[["Boot"]]), info = m)
    if (startsWith(m, "FE-")) expect_true(is.finite(p[["Asy"]]), info = m)
    else expect_false("Asy" %in% names(p), info = m)
  }
})

test_that("cluster-level and within-cluster exposures need no shuffle mode", {
  d <- make_input()
  ids <- lapply(d$X, function(x)
    factor(rep(1:20, each = 2), levels = 1:25))
  between <- MetaQCAT(d$counts, d$X, 1, cluster.id = ids,
                      taxonomy = d$taxonomy, method = "FE-MetaQCAT-O",
                      n.boot = 9L)
  expect_length(between$failed.lineage, 0L)
  expect_true(all(is.finite(between$lineage.pval)))
  expect_equal(unname(between$bootstrap$n.clusters), c(20L, 20L))

  within.x <- lapply(d$X, function(x) {
    x[, 1] <- rep(0:1, 20)
    x
  })
  within <- MetaQCAT(d$counts, within.x, 1, cluster.id = ids,
                     taxonomy = d$taxonomy, method = "RE-MetaQCAT-O",
                     n.boot = 9L)
  expect_length(within$failed.lineage, 0L)
  expect_true(all(is.finite(within$lineage.pval)))
})

test_that("study and taxon alignment preserve lineage dimensions", {
  d <- make_input()
  ref <- MetaQCAT(d$counts, d$X, 1, taxonomy = d$taxonomy,
                  n.boot = 9L, seed = 17L)
  y <- d$counts
  y[[2]] <- y[[2]][, 6:1]
  fit <- MetaQCAT(y, rev(d$X), 1, taxonomy = d$taxonomy[6:1, ],
                  n.boot = 9L, seed = 17L)
  expect_equal(fit$lineage.pval, ref$lineage.pval)
  tax <- d$taxonomy[1:3, , drop = FALSE]
  y <- lapply(d$counts, function(a) a[, 1:3])
  one <- MetaQCAT(y, d$X, 1, taxonomy = tax, n.boot = 9L)
  expect_identical(dim(one$lineage.pval), c(2L, 1L))
  expect_identical(colnames(one$lineage.pval), "FamilyA")
})

test_that("failed lineages retain their NA position and reason", {
  d <- make_input()
  y <- lapply(d$counts, function(a) {a[, 1:3] <- 0; a})
  expect_warning(
    fit <- MetaQCAT(y, d$X, 1, taxonomy = d$taxonomy,
                    method = "FE-MetaQCAT", n.boot = 9L),
    "NA columns retained"
  )
  expect_identical(colnames(fit$lineage.pval), c("FamilyA", "FamilyB"))
  expect_true(all(is.na(fit$lineage.pval[, "FamilyA"])))
  expect_true(all(is.finite(fit$lineage.pval[, "FamilyB"])))
  expect_identical(names(fit$failed.lineage), "FamilyA")
  expect_false("FamilyA" %in% fit$sig.lineage)
})

test_that("depth filtering excludes only empty studies", {
  d <- make_input()
  y <- d$counts
  y[[1]][seq(1L, 39L, by = 2L), ] <- 0
  fit <- MetaQCAT(y, d$X, 1, n.boot = 9L)
  expect_null(fit$error)
  y[[1]][, ] <- 0
  expect_warning(fit <- MetaQCAT(y, d$X, 1, n.boot = 9L),
                 "Excluded empty studies")
  expect_null(fit$error)
  y[[2]][, ] <- 0
  expect_error(MetaQCAT(y, d$X, 1, n.boot = 9L), "No studies remain")
})

test_that("cluster bootstrap is reproducible across worker schedules", {
  skip_on_os("windows")
  d <- make_input()
  ids <- lapply(d$X, function(x) rep(seq_len(20), each = 2))
  a <- MetaQCAT(d$counts, d$X, 1, taxonomy = d$taxonomy,
                cluster.id = ids, n.boot = 19L, seed = 7L,
                method = "RE-MetaQCAT-O")
  b <- MetaQCAT(d$counts, d$X, 1, taxonomy = d$taxonomy,
                cluster.id = ids, n.boot = 19L, seed = 7L,
                method = "RE-MetaQCAT-O", n.threads = 2L)
  expect_identical(a$lineage.pval, b$lineage.pval)
})

test_that("multiple tested covariates retain coefficient dimensions", {
  d <- make_input()
  x <- lapply(d$X, function(a)
    cbind(a, trend = seq_len(nrow(a)) %% 3))
  fit <- MetaQCAT(d$counts, x, 1:2, taxonomy = d$taxonomy,
                  method = "FE-MetaQCAT", n.boot = 9L,
                  save.coef = TRUE)
  expect_length(fit$failed.lineage, 0L)
  expect_true(all(is.finite(fit$lineage.pval)))
  r <- fit$lineage.coef[[1L]]
  if (!is.null(r$est.beta.comp)) expect_identical(nrow(r$est.beta.comp), 2L)
})

test_that("packaged CRC and IBD designs remain distinct", {
  data("CRC_data", package = "MetaQCAT", envir = environment())
  data("IBD_data", package = "MetaQCAT", envir = environment())
  expect_false("design" %in% names(CRC_data))
  expect_null(CRC_data$subject_ids)
  expect_true("design" %in% names(IBD_data))
  expect_identical(sum(vapply(CRC_data$counts, nrow, integer(1))), 573L)
  expect_identical(sum(vapply(IBD_data$counts, nrow, integer(1))), 3141L)
  members <- do.call(rbind, lapply(names(IBD_data$metadata), function(study) {
    m <- IBD_data$metadata[[study]]
    unique(data.frame(subject = paste(m$dataset_name, m$subject_accession,
                                      sep = "::"), study = study))
  }))
  expect_false(anyDuplicated(members$subject) > 0L)
})
