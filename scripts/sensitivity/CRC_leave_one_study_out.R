# Run from the MetaQCAT-repro repository root.
source(file.path("R", "project_paths.R"))

# Leave-one-study-out sensitivity analysis for the CRC bootstrap MetaQCAT-O runs.
# Writes only to res/CRC/wild_score_bootstrap/leave_one_study_out.
script.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
# Repository paths are defined by R/project_paths.R.

library(MetaQCAT)
stopifnot(packageVersion("MetaQCAT") == "0.2.1")

data.dir <- file.path(script.dir, "data", "CRC")
result.dir <- file.path(script.dir, "res", "CRC", "wild_score_bootstrap",
                        "leave_one_study_out")
dir.create(result.dir, recursive = TRUE, showWarnings = FALSE)
n.boot <- 1000L
n.threads <- 7L
workers.used <- min(n.threads, parallelly::availableCores())
seed <- 1234L
min.depth <- 0
# Match the default used by CRC_all_methods.R.
min.prev <- 0.1
max.zero.fraction <- 0.90
depth.threshold <- 1000
methods <- c("FE-MetaQCAT-O", "RE-MetaQCAT-O")

# Reuse the exact prepared study tables used by CRC_all_methods.R.
# Study omission does not redo the global depth or genus screening.
processed.env <- new.env(parent = emptyenv())
load(file.path(script.dir, "res", "CRC", "wild_score_bootstrap",
               "CRC_processed.RData"), envir = processed.env)
counts <- processed.env$counts
X <- processed.env$disease
taxonomy <- processed.env$taxonomy[, c("Rank5", "Rank6"), drop = FALSE]
colnames(taxonomy) <- c("Rank1", "Rank2")
stopifnot(identical(names(counts), names(X)),
          length(counts) == 5L,
          identical(processed.env$preprocess.settings$min.depth, 1000L),
          identical(processed.env$preprocess.settings$max.zero.fraction, 0.90),
          identical(processed.env$preprocess.settings$min.prev, min.prev))

primary.env <- new.env(parent = emptyenv())
load(file.path(script.dir, "res", "CRC", "wild_score_bootstrap",
               "CRC_MetaQCAT_bootstrap_res.RData"),
     envir = primary.env)
primary.results <- list("FE-MetaQCAT-O" = primary.env$FE_CRC,
                        "RE-MetaQCAT-O" = primary.env$RE_CRC)
expected.lineages <- colnames(primary.results[[1L]]$lineage.qval)

extract.lineages <- function(result, method, omitted.study) {
  p.row <- paste0("Bootstrap-", method)
  p <- q <- setNames(rep(NA_real_, length(expected.lineages)), expected.lineages)
  reason <- setNames(rep("", length(expected.lineages)), expected.lineages)
  if (inherits(result, "error")) {
    reason[] <- conditionMessage(result)
  } else if (!is.null(result$lineage.pval) && p.row %in% rownames(result$lineage.pval)) {
    found <- intersect(expected.lineages, colnames(result$lineage.pval))
    p[found] <- result$lineage.pval[p.row, found]
    q[found] <- result$lineage.qval[p.row, found]
    failed <- result$failed.lineage
    if (length(failed)) {
      common <- intersect(names(failed), expected.lineages)
      reason[common] <- unname(failed[common])
    }
    reason[setdiff(expected.lineages, colnames(result$lineage.pval))] <-
      "Lineage not returned by this fit."
  } else {
    reason[] <- "Requested bootstrap p-value row was not returned."
  }
  data.frame(omitted_study = omitted.study, method = method,
    lineage = expected.lineages, p = unname(p), q = unname(q),
    significant = ifelse(is.na(q), NA, q < 0.05),
    failure_reason = unname(reason), stringsAsFactors = FALSE)
}

lineage.rows <- list()
run.rows <- list()
for (method in methods)
  lineage.rows[[paste("full", method, sep = "_")]] <-
    extract.lineages(primary.results[[method]], method, "None")

for (study in names(counts)) {
  keep.study <- setdiff(names(counts), study)
  n.samples <- sum(vapply(counts[keep.study], nrow, integer(1)))
  n.cases <- sum(vapply(keep.study, function(s) sum(X[[s]][, 1L] == 1), numeric(1)))
  n.controls <- sum(vapply(keep.study, function(s) sum(X[[s]][, 1L] == 0), numeric(1)))
  message("CRC LOSO: omitting ", study, "; ", n.samples, " samples remain.")
  for (method in methods) {
    started <- Sys.time()
    fit <- tryCatch(MetaQCAT::MetaQCAT(
      counts = counts[keep.study], X = X[keep.study], X.index = 1L,
      taxonomy = taxonomy, method = method, min.depth = min.depth,
      min.prev = min.prev, n.boot = n.boot, n.threads = n.threads,
      seed = seed), error = function(e) e)
    seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    safe.study <- gsub("[^A-Za-z0-9_.-]", "_", study)
    saveRDS(fit, file.path(result.dir,
      paste0("omit_", safe.study, "_", gsub("[^A-Za-z0-9]+", "_", method), ".rds")))
    lineage.rows[[paste(study, method, sep = "_")]] <-
      extract.lineages(fit, method, study)
    run.rows[[paste(study, method, sep = "_")]] <- data.frame(
      omitted_study = study, method = method, remaining_studies = length(keep.study),
      samples = n.samples, cases = n.cases, controls = n.controls,
      bootstrap_budget = n.boot, seed = seed, requested_threads = n.threads,
      workers_used = workers.used, seconds = seconds,
      failed_lineages = if (inherits(fit, "error")) NA_integer_ else
        length(fit$failed.lineage),
      error = if (inherits(fit, "error")) conditionMessage(fit) else "",
      stringsAsFactors = FALSE)
    saveRDS(list(lineages = do.call(rbind, lineage.rows),
      runs = do.call(rbind, run.rows)), file.path(result.dir, "CRC_loso_checkpoint.rds"))
    utils::write.csv(do.call(rbind, lineage.rows),
      file.path(result.dir, "CRC_loso_lineage_results.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, run.rows),
      file.path(result.dir, "CRC_loso_run_summary.csv"), row.names = FALSE)
  }
}

all.lineages <- do.call(rbind, lineage.rows)
full.lineages <- all.lineages[all.lineages$omitted_study == "None", , drop = FALSE]
stability.rows <- list()
for (study in c("None", names(counts))) for (method in methods) {
  current <- all.lineages[all.lineages$omitted_study == study &
    all.lineages$method == method, , drop = FALSE]
  baseline <- full.lineages[full.lineages$method == method, , drop = FALSE]
  primary.sig <- baseline$lineage[!is.na(baseline$q) & baseline$q < 0.05]
  current.sig <- current$lineage[!is.na(current$q) & current$q < 0.05]
  retained <- intersect(primary.sig, current.sig)
  union.sig <- union(primary.sig, current.sig)
  stability.rows[[paste(study, method, sep = "_")]] <- data.frame(
    omitted_study = study, method = method,
    primary_discoveries = length(primary.sig), discoveries = length(current.sig),
    retained = length(retained), lost = length(setdiff(primary.sig, current.sig)),
    gained = length(setdiff(current.sig, primary.sig)),
    retention_fraction = if (!length(primary.sig)) NA_real_ else
      length(retained) / length(primary.sig),
    jaccard = if (!length(union.sig)) NA_real_ else
      length(retained) / length(union.sig),
    retained_lineages = paste(retained, collapse = ";"),
    lost_lineages = paste(setdiff(primary.sig, current.sig), collapse = ";"),
    gained_lineages = paste(setdiff(current.sig, primary.sig), collapse = ";"),
    stringsAsFactors = FALSE)
}

utils::write.csv(all.lineages,
  file.path(result.dir, "CRC_loso_lineage_results.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, run.rows),
  file.path(result.dir, "CRC_loso_run_summary.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, stability.rows),
  file.path(result.dir, "CRC_loso_stability_summary.csv"), row.names = FALSE)
saveRDS(list(settings = list(n.boot = n.boot, n.threads = n.threads,
  seed = seed, min.depth = min.depth, min.prev = min.prev,
  requested_threads = n.threads, workers_used = workers.used,
  max.zero.fraction = max.zero.fraction, depth.threshold = depth.threshold),
  lineage_results = all.lineages, run_summary = do.call(rbind, run.rows),
  stability_summary = do.call(rbind, stability.rows)),
  file.path(result.dir, "CRC_loso_results.rds"))
message("CRC LOSO complete: ", result.dir)
