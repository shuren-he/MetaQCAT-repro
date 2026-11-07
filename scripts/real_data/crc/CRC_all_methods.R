# CRC: one shared preprocessing step, MetaQCAT, and all P/CD/RD comparators.
# Run from the repository root, or execute the numbered sections in RStudio.
# This is a fresh analysis. Running it replaces the corresponding outputs in
# this repository; it does not reuse retained historical comparator values.
source(file.path("R", "project_paths.R"))

## 1. Configuration -----------------------------------------------------------
min.depth <- 1000L
max.zero.fraction <- 0.90
min.prev <- 0.10                 # MetaQCAT component-specific prevalence filter
n.boot <- 1000L                  # MetaQCAT wild score-bootstrap draws
n.perm <- 1000L                  # Comparator label permutations
n.threads <- 7L                  # MetaQCAT workers; comparator loop is explicit
base.seed <- 1234L
fdr.alpha <- 0.05

for (package in c("MetaQCAT", "miLineage", "HMP", "MMUPHin")) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install the analysis dependencies first: missing ", package)
}
stopifnot(packageVersion("MetaQCAT") == "0.2.1")
data.dir <- file.path(project.dir, "data", "CRC")
result.dir <- file.path(project.dir, "res", "CRC")
bootstrap.dir <- file.path(result.dir, "wild_score_bootstrap")
table.dir <- file.path(project.dir, "tables")

## 2. Load once and align sample identifiers ----------------------------------
count.env <- new.env(parent = emptyenv())
meta.env <- new.env(parent = emptyenv())
tax.env <- new.env(parent = emptyenv())
load(file.path(data.dir, "count.genus.RData"), envir = count.env)
load(file.path(data.dir, "meta.RData"), envir = meta.env)
load(file.path(data.dir, "tax.RData"), envir = tax.env)
count.matrix <- as.matrix(count.env$count.genus)
sample.metadata <- as.data.frame(meta.env$meta, stringsAsFactors = FALSE)
source.taxonomy <- tax.env$tax
stopifnot(nrow(count.matrix) == nrow(sample.metadata),
          !anyNA(count.matrix), all(is.finite(count.matrix)),
          all(count.matrix >= 0 & count.matrix == floor(count.matrix)),
          !anyDuplicated(rownames(count.matrix)),
          !anyDuplicated(colnames(count.matrix)),
          !anyNA(sample.metadata$Sample_ID),
          !anyDuplicated(sample.metadata$Sample_ID))
count.keys <- make.names(sample.metadata$Sample_ID)
stopifnot(!anyDuplicated(count.keys),
          setequal(rownames(count.matrix), count.keys))
sample.metadata <- sample.metadata[
  match(rownames(count.matrix), count.keys), , drop = FALSE]
# Preserve the original count-table analysis keys and row order. The original
# Sample_ID is still available in metadata for exclusion reports.
rownames(sample.metadata) <- rownames(count.matrix)

## 3. ONE depth filter shared by EVERY method ---------------------------------
keep <- rowSums(count.matrix) >= min.depth
excluded.samples <- data.frame(
  sample_id = sample.metadata$Sample_ID[!keep],
  study = sample.metadata$Country[!keep],
  group = sample.metadata$Group[!keep],
  original_depth = rowSums(count.matrix)[!keep],
  stringsAsFactors = FALSE)
count.matrix <- count.matrix[keep, , drop = FALSE]
sample.metadata <- sample.metadata[keep, , drop = FALSE]
stopifnot(identical(rownames(count.matrix), rownames(sample.metadata)),
          !anyNA(sample.metadata$Country),
          all(sample.metadata$Group %in% c("CRC", "CTR")))

## 4. Taxonomy, study tables, and common settings -------------------------------
ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
taxonomy <- source.taxonomy[, ranks, drop = FALSE]
taxonomy <- as.matrix(taxonomy[!duplicated(taxonomy$genus), , drop = FALSE])
rownames(taxonomy) <- taxonomy[, "genus"]
colnames(taxonomy) <- paste0("Rank", seq_len(ncol(taxonomy)))
stopifnot(all(colnames(count.matrix) %in% rownames(taxonomy)))
family.taxonomy <- taxonomy[, c("Rank5", "Rank6"), drop = FALSE]
colnames(family.taxonomy) <- c("Rank1", "Rank2")

metadata <- split(sample.metadata, sample.metadata$Country, drop = TRUE)
counts <- lapply(metadata, function(m) {
  y <- count.matrix[rownames(m), , drop = FALSE]
  y[, colMeans(y == 0) <= max.zero.fraction, drop = FALSE]
})
X <- lapply(metadata, function(m) {
  matrix(as.integer(m$Group == "CRC"), ncol = 1L,
         dimnames = list(rownames(m), "disease"))
})
disease <- X
stopifnot(all(vapply(names(counts), function(study)
  identical(rownames(counts[[study]]), rownames(X[[study]])), logical(1))),
  all(vapply(counts, function(y) ncol(y) > 0L && all(rowSums(y) > 0),
             logical(1))))
study.summary <- data.frame(
  study = names(counts),
  samples = vapply(counts, nrow, integer(1)),
  cases = vapply(X, function(x) sum(x[, 1L] == 1L), integer(1)),
  controls = vapply(X, function(x) sum(x[, 1L] == 0L), integer(1)),
  retained_genera = vapply(counts, ncol, integer(1)), row.names = NULL)
print(study.summary)

# CD/RD pool the SAME retained samples. Their pooled taxon filter uses the
# same threshold, applied to the pooled original matrix (not zero-filled
# versions of the study-filtered tables). This preserves the P/CD/RD strategy.
pooled.metadata <- do.call(rbind, unname(metadata))
pooled <- count.matrix[rownames(pooled.metadata), , drop = FALSE]
pooled <- pooled[, colMeans(pooled == 0) <= max.zero.fraction, drop = FALSE]
pooled.group <- setNames(as.integer(pooled.metadata$Group == "CRC"),
                         rownames(pooled.metadata))
stopifnot(setequal(rownames(pooled), unlist(lapply(counts, rownames))),
          identical(rownames(pooled), rownames(pooled.metadata)))

## 5. MetaQCAT: two fits expose all six component/omnibus tests ------------------
dir.create(bootstrap.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table.dir, recursive = TRUE, showWarnings = FALSE)
subject_ids <- NULL
adjust <- lapply(X, function(x)
  matrix(numeric(), nrow(x), 0L, dimnames = list(rownames(x), NULL)))
preprocess.settings <- list(min.depth = min.depth,
  max.zero.fraction = max.zero.fraction, min.prev = min.prev, seed = base.seed)
save(counts, metadata, disease, adjust, subject_ids, taxonomy,
     preprocess.settings,
     file = file.path(bootstrap.dir, "CRC_processed.RData"), compress = "xz")
run.started <- Sys.time()
meta.started <- Sys.time()
FE_CRC <- MetaQCAT::MetaQCAT(
  counts = counts, X = X, X.index = 1L, taxonomy = family.taxonomy,
  method = "FE-MetaQCAT-O", min.depth = 0, min.prev = min.prev,
  n.boot = n.boot, n.threads = n.threads, seed = base.seed)
fe.seconds <- as.numeric(difftime(Sys.time(), meta.started, units = "secs"))
meta.started <- Sys.time()
RE_CRC <- MetaQCAT::MetaQCAT(
  counts = counts, X = X, X.index = 1L, taxonomy = family.taxonomy,
  method = "RE-MetaQCAT-O", min.depth = 0, min.prev = min.prev,
  n.boot = n.boot, n.threads = n.threads, seed = base.seed)
re.seconds <- as.numeric(difftime(Sys.time(), meta.started, units = "secs"))
save(FE_CRC, RE_CRC, file = file.path(bootstrap.dir,
     "CRC_MetaQCAT_bootstrap_res.RData"), compress = "xz")

# Derive the family universe from the new fit, not an old result-table template.
families <- colnames(FE_CRC$lineage.pval)
stopifnot(setequal(families, colnames(RE_CRC$lineage.pval)))
family.keys <- paste0(families, " (family)")
our.columns <- c("FE.MetaQCAT.asy", "FE.MetaQCAT*.asy", "FE.MetaQCAT.O.asy",
                 "FE.MetaQCAT.res", "FE.MetaQCAT*.res", "FE.MetaQCAT.O.res",
                 "RE.MetaQCAT", "RE.MetaQCAT*", "RE.MetaQCAT.O")
fields <- rep(c("pvalue.comp", "pvalue.pa", "pvalue.comb"), 3L)
calibrations <- c(rep("Asy", 3L), rep("Boot", 6L))
our.p <- matrix(NA_real_, length(families), length(our.columns),
                dimnames = list(family.keys, our.columns))
for (j in seq_along(our.columns)) {
  fit <- if (j <= 6L) FE_CRC else RE_CRC
  our.p[, j] <- vapply(fit$lineage.coef[families], function(z) {
    p <- z[[fields[j]]]
    if (is.null(p) || !calibrations[j] %in% names(p)) NA_real_ else
      unname(p[calibrations[j]])
  }, numeric(1))
}

## 6. Comparator helpers (existing statistical calculations retained) -----------
rarefy.counts <- function(y, seed) {
  y <- as.matrix(y)
  if (anyNA(y) || any(y < 0) || any(y != floor(y))) {
    stop("Rarefaction requires nonnegative integer counts.")
  }
  target <- as.integer(min(rowSums(y)))
  if (target < 1L) stop("A sample has no reads after filtering.")
  set.seed(seed)
  result <- t(vapply(seq_len(nrow(y)), function(i) {
    tabulate(sample(rep.int(seq_len(ncol(y)), as.integer(y[i, ])),
                    size = target, replace = FALSE), nbins = ncol(y))
  }, integer(ncol(y))))
  dimnames(result) <- dimnames(y)
  result
}

fisher.p <- function(p) {
  p <- as.numeric(p)
  p <- p[is.finite(p) & p >= 0 & p <= 1]
  if (length(p) < 2L) return(NA_real_)
  if (any(p == 0)) return(0)
  stats::pchisq(-2 * sum(log(p)),
                df = 2 * length(p), lower.tail = FALSE)
}

run.method <- function(y, y.rare, group, strategy, method, seed) {
  y <- if (method == "QCAT2") y.rare else y
  y <- y[, colSums(y) > 0, drop = FALSE]
  # Only historical P-QCAT1 used min.depth=1 before its Tax=NULL call.
  # For the other miLineage methods, use a one-family Tax matrix and let the
  # installed package's lineage path handle family-zero samples itself.
  if (strategy == "P" && method == "QCAT1") {
    y <- y[rowSums(y) > 0, , drop = FALSE]
  }
  if (ncol(y) < 2L) stop("Fewer than two observed genera.")
  group <- as.integer(group[rownames(y)])
  if (anyNA(group) || any(tabulate(group + 1L, nbins = 2L) < 1L)) {
    stop("Case/control labels are missing or one group is empty.")
  }
  x <- matrix(group, ncol = 1L, dimnames = list(rownames(y), "disease"))
  family.name <- unique(taxonomy[colnames(y), "Rank5"])
  if (length(family.name) != 1L) stop("Taxonomy does not identify one family.")
  family.tax <- cbind(Rank1 = rep(family.name, ncol(y)),
                      Rank2 = colnames(y))
  rownames(family.tax) <- colnames(y)
  set.seed(seed)
  if (method == "QCAT1") {
    if (strategy == "P") {
      answer <- miLineage::QCAT(y, x, X.index = 1L, Tax = NULL,
                               min.depth = 0, n.perm = n.perm)$pval
    } else {
      answer <- miLineage::QCAT(y, x, X.index = 1L, Tax = family.tax,
                               min.depth = 0, n.perm = n.perm)$lineage.pval[, 1L]
    }
    n.success <- NA_integer_  # Not exposed by miLineage.
  } else if (method == "QCAT2") {
    answer <- miLineage::QCAT_GEE(y, x, X.index = 1L, Z = x, Z.index = 1L,
                                 Tax = family.tax, min.depth = 0,
                                 n.perm = n.perm)$lineage.pval[["Two-Part"]][, 1L]
    n.success <- NA_integer_
  } else if (method == "ZIGDM") {
    answer <- miLineage::ZIGDM(y, X4mean = x, X4disp = NULL, X4zero = NULL,
                              test.type = "Mean", X.index = 1L,
                              Tax = family.tax, min.depth = 0,
                              n.perm = n.perm)$lineage.pval[, 1L]
    n.success <- NA_integer_
  } else if (method == "DM") {
    grouped <- split.data.frame(y, factor(group, levels = 0:1))
    names(grouped) <- c("control", "case")
    observed <- HMP::Xmcupo.sevsample(grouped)
    statistic <- as.numeric(observed$`Xmcupo statistics`)
    asymptotic <- as.numeric(observed$`p value`)
    if (length(statistic) != 1L || !is.finite(statistic)) {
      stop("Observed DM statistic is unavailable.")
    }
    exceed <- 0L
    n.success <- 0L
    for (b in seq_len(n.perm)) {
      permuted.group <- sample(group)
      groups <- split.data.frame(y, factor(permuted.group, levels = 0:1))
      names(groups) <- c("control", "case")
      permuted <- tryCatch(HMP::Xmcupo.sevsample(groups)$`Xmcupo statistics`,
                           error = function(e) NA_real_)
      if (length(permuted) == 1L && is.finite(permuted)) {
        n.success <- n.success + 1L
        exceed <- exceed + as.integer(permuted >= statistic)
      }
    }
    # This DM comparator uses exceed / n.success (no plus-one). Never
    # present a partial run as a 1,000-permutation result.
    resampling <- if (n.success == n.perm) {
      exceed / n.success
    } else {
      NA_real_
    }
    answer <- c(asymptotic, resampling)
  } else {
    stop("Unknown comparator method: ", method)
  }
  answer <- as.numeric(answer)
  if (length(answer) != 2L) {
    stop("Method did not return two p-values.")
  }
  answer[!is.finite(answer) | answer < 0 | answer > 1] <- NA_real_
  if (all(is.na(answer))) stop("Method returned no valid p-values.")
  list(p = answer, n_perm_success = n.success,
       n_samples = if (method == "ZIGDM") sum(rowSums(y) > 0) else nrow(y),
       n_taxa = ncol(y))
}

record.test <- function(family, strategy, study, method, y, y.rare,
                        group, seed) {
  started <- proc.time()[["elapsed"]]
  warnings <- character()
  result <- withCallingHandlers(
    tryCatch(run.method(y, y.rare, group, strategy, method, seed),
             error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  failed <- inherits(result, "error")
  detail <- if (failed) conditionMessage(result) else if (anyNA(result$p)) {
    paste0("Partial p-values; ", result$n_perm_success,
           " of ", n.perm, " DM permutations succeeded or method returned NA.")
  } else ""
  if (length(warnings)) {
    detail <- paste(c(detail, unique(warnings)[seq_len(min(3L, length(unique(warnings))))]),
                    collapse = " | ")
  }
  row <- data.frame(
    family = family, strategy = strategy, study = study, method = method,
    n_samples = if (failed) NA_integer_ else result$n_samples,
    n_taxa = if (failed) NA_integer_ else result$n_taxa,
    n_perm_requested = n.perm,
    n_perm_success = if (failed) NA_integer_ else result$n_perm_success,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    status = if (failed) "failed" else if (anyNA(result$p)) "partial" else "ok",
    detail = detail,
    stringsAsFactors = FALSE
  )
  list(p = if (failed) c(NA_real_, NA_real_) else result$p,
       diagnostic = row)
}


## 7. P/CD/RD: batch correction and full-composition rarefaction ONCE ------------
methods <- c("QCAT1", "QCAT2", "DM", "ZIGDM")
strategies <- c("P", "CD", "RD")
column.name <- function(strategy, method, component) {
  prefix <- if (strategy == "P") "PC" else strategy
  method.label <- if (strategy != "P" && method == "QCAT1") "QCAT" else method
  paste(prefix, method.label, component, sep = ".")
}
baseline.columns <- unlist(lapply(strategies, function(strategy) {
  unlist(lapply(c("asy", "res"), function(component)
    vapply(methods, function(method) column.name(strategy, method, component), "")))
}), use.names = FALSE)
baseline.p <- matrix(NA_real_, length(families), length(baseline.columns),
                     dimnames = list(family.keys, baseline.columns))
diagnostics <- data.frame(
  family = character(), strategy = character(), study = character(),
  method = character(), n_samples = integer(), n_taxa = integer(),
  n_perm_requested = integer(), n_perm_success = integer(),
  elapsed_seconds = numeric(), status = character(), detail = character(),
  stringsAsFactors = FALSE)

batch.fit <- MMUPHin::adjust_batch(
  feature_abd = t(pooled), batch = "Country", covariates = "Group",
  data = pooled.metadata, control = list(verbose = FALSE))
corrected <- t(batch.fit$feature_abd_adj)
corrected <- corrected[rownames(pooled), colnames(pooled), drop = FALSE]
stopifnot(all(is.finite(corrected)), all(corrected >= 0))
rarefied <- list(
  P = lapply(seq_along(counts), function(i)
    rarefy.counts(counts[[i]], base.seed + i)),
  CD = rarefy.counts(corrected, base.seed + 100L),
  RD = rarefy.counts(pooled, base.seed + 200L))
names(rarefied$P) <- names(counts)

## 8. Explicit family/method loops (RStudio-friendly, no hidden runner) ----------
for (family.index in seq_along(families)) {
  family <- families[family.index]
  family.key <- family.keys[family.index]
  family.taxa <- rownames(taxonomy)[taxonomy[, "Rank5"] == family]
  message("CRC family ", family.index, "/", length(families), ": ", family)

  for (method.index in seq_along(methods)) {
    method <- methods[method.index]
    p.by.study <- matrix(NA_real_, length(counts), 2L,
      dimnames = list(names(counts), c("asy", "res")))
    for (study.index in seq_along(counts)) {
      study <- names(counts)[study.index]
      taxa <- intersect(colnames(counts[[study]]), family.taxa)
      y <- counts[[study]][, taxa, drop = FALSE]
      y.rare <- rarefied$P[[study]][, taxa, drop = FALSE]
      group <- setNames(as.integer(X[[study]][, 1L]), rownames(X[[study]]))
      answer <- record.test(family, "P", study, method, y, y.rare, group,
        base.seed + 10000L * family.index + 100L * method.index + study.index)
      p.by.study[study.index, ] <- answer$p
      diagnostics <- rbind(diagnostics, answer$diagnostic)
    }
    baseline.p[family.key, column.name("P", method, "asy")] <-
      fisher.p(p.by.study[, 1L])
    baseline.p[family.key, column.name("P", method, "res")] <-
      fisher.p(p.by.study[, 2L])
  }

  for (strategy.index in seq_along(c("CD", "RD"))) {
    strategy <- c("CD", "RD")[strategy.index]
    y.full <- if (strategy == "CD") corrected else pooled
    taxa <- intersect(colnames(y.full), family.taxa)
    y <- y.full[, taxa, drop = FALSE]
    y.rare <- rarefied[[strategy]][, taxa, drop = FALSE]
    for (method.index in seq_along(methods)) {
      method <- methods[method.index]
      answer <- record.test(family, strategy, "pooled", method, y, y.rare,
        pooled.group, base.seed + 10000L * family.index +
          1000L * strategy.index + 100L * method.index)
      baseline.p[family.key, c(column.name(strategy, method, "asy"),
                              column.name(strategy, method, "res"))] <- answer$p
      diagnostics <- rbind(diagnostics, answer$diagnostic)
    }
  }
  saveRDS(list(p_values = baseline.p, diagnostics = diagnostics,
               completed_families = families[seq_len(family.index)]),
          file.path(result.dir, "CRC_all_methods_partial.rds"), compress = "xz")
}

## 9. Fresh p/q tables; no pre-existing comparator result/template required -----
p.values <- cbind(our.p, baseline.p)
q.values <- apply(p.values, 2L, p.adjust, method = "BH")
resampling.columns <- c(our.columns[4:9],
                        baseline.columns[endsWith(baseline.columns, ".res")])
write.csv(p.values, file.path(result.dir, "family_p_val.csv"))
write.csv(q.values, file.path(result.dir, "family_q_val.csv"))
write.csv(q.values[, resampling.columns, drop = FALSE],
          file.path(result.dir, "family_q_val_resampling.csv"))
paper.q <- q.values
colnames(paper.q)[1:9] <- c(
  "FE-MetaQCAT (asy)", "FE-MetaQCAT* (asy)", "FE-MetaQCAT-O (asy)",
  "FE-MetaQCAT (boot)", "FE-MetaQCAT* (boot)", "FE-MetaQCAT-O (boot)",
  "RE-MetaQCAT (boot)", "RE-MetaQCAT* (boot)", "RE-MetaQCAT-O (boot)")
write.csv(paper.q, file.path(table.dir, "CRC_family_q_val.csv"))

# The comparison UpSet uses resampling q-values only: a union of the four
# comparator methods per P/CD/RD strategy, and the FE/RE omnibus tests.
sig <- data.frame(Lineage = family.keys, check.names = FALSE)
for (strategy in c("RD", "P", "CD")) {
  prefix <- if (strategy == "P") "PC." else paste0(strategy, ".")
  cols <- baseline.columns[startsWith(baseline.columns, prefix) &
                             endsWith(baseline.columns, ".res")]
  sig[[strategy]] <- as.integer(rowSums(
    q.values[, cols, drop = FALSE] < fdr.alpha, na.rm = TRUE) > 0)
}
sig$FE <- as.integer(!is.na(q.values[, "FE.MetaQCAT.O.res"]) &
                       q.values[, "FE.MetaQCAT.O.res"] < fdr.alpha)
sig$RE <- as.integer(!is.na(q.values[, "RE.MetaQCAT.O"]) &
                       q.values[, "RE.MetaQCAT.O"] < fdr.alpha)
write.csv(sig, file.path(result.dir, "sig_family.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(result.dir, "CRC_all_methods_diagnostics.csv"),
          row.names = FALSE)
write.csv(study.summary, file.path(result.dir, "CRC_all_methods_study_summary.csv"),
          row.names = FALSE)
write.csv(excluded.samples, file.path(result.dir, "CRC_all_methods_exclusions.csv"),
          row.names = FALSE)
stopifnot(isTRUE(all.equal(unname(q.values[, "FE.MetaQCAT.O.res"]),
  unname(FE_CRC$lineage.qval["Bootstrap-FE-MetaQCAT-O", families]))),
  isTRUE(all.equal(unname(q.values[, "RE.MetaQCAT.O"]),
  unname(RE_CRC$lineage.qval["Bootstrap-RE-MetaQCAT-O", families]))))

## 10. Save the complete analysis object ----------------------------------------
CRC_all_methods <- list(
  FE_CRC = FE_CRC, RE_CRC = RE_CRC,
  p_values = p.values, q_values = q.values, significant = sig,
  diagnostics = diagnostics, study_summary = study.summary,
  excluded_samples = excluded.samples,
  settings = list(min.depth = min.depth, max.zero.fraction = max.zero.fraction,
    min.prev = min.prev, n.boot = n.boot, n.perm = n.perm,
    seed = base.seed, n.threads = n.threads,
    pooled_filter = "same threshold applied to pooled original counts",
    dm_permutation_p = "exceed/success; all requested draws required",
    input_md5 = tools::md5sum(file.path(data.dir,
      c("count.genus.RData", "meta.RData", "tax.RData"))),
    script_md5 = unname(tools::md5sum(file.path(project.dir, "scripts",
      "real_data", "crc", "CRC_all_methods.R"))),
    package.versions = vapply(c("MetaQCAT", "miLineage", "HMP", "MMUPHin"),
      function(p) as.character(packageVersion(p)), character(1))),
  elapsed_seconds = as.numeric(difftime(Sys.time(), run.started, units = "secs")),
  meta_seconds = c(FE = fe.seconds, RE = re.seconds),
  session_info = sessionInfo())
saveRDS(CRC_all_methods, file.path(result.dir, "CRC_all_methods.rds"),
        compress = "xz")
message("Saved CRC analyses and fresh result tables under ", result.dir)
