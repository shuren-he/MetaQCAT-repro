# Run from the MetaQCAT-repro repository root.
source(file.path("R", "project_paths.R"))

# IBD leave-one-source-cohort-out analysis using the saved primary scores.
# Rarefaction, study-specific eligibility, and reference taxa remain fixed.
# No null model is refitted. All strata from a source cohort are omitted together.
script.started <- Sys.time()
script.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
# Repository paths are defined by R/project_paths.R.
library(MetaQCAT)
stopifnot(as.character(utils::packageVersion("MetaQCAT")) == "0.2.1")

primary.dir <- file.path(script.dir, "res", "IBD", "wild_score_bootstrap")
primary.file <- file.path(primary.dir, "IBD_primary.rds")
data.file <- file.path(primary.dir, "IBD_processed_bootstrap.RData")
result.dir <- file.path(primary.dir, "leave_one_cohort_out")
dir.create(result.dir, recursive = TRUE, showWarnings = FALSE)
n.boot <- 300L
seed <- 1234L
methods <- c("FE-MetaQCAT-O", "RE-MetaQCAT-O")
data.env <- new.env(parent = emptyenv())
load(data.file, envir = data.env)
counts <- data.env$counts
design <- data.env$design
subject_ids <- data.env$subject_ids
metadata <- data.env$metadata
primary <- readRDS(primary.file)$results[methods]
lineages <- colnames(primary[[1L]]$lineage.pval)
source.md5 <- unname(tools::md5sum(primary.file))
source.cohort <- vapply(metadata, function(m) {
  value <- unique(as.character(m$dataset_name))
  if (length(value) != 1L) stop("A stratum maps to multiple source cohorts.")
  value
}, character(1))
cohorts <- unique(unname(source.cohort))
stopifnot(length(lineages) == 45L, length(cohorts) == 7L,
          identical(names(counts), names(design)),
          identical(names(counts), names(subject_ids)),
          identical(names(counts), names(metadata)))

# FE and RE fitted the same reduced models. Verify this before reusing them.
for (lineage in lineages) {
  fe <- primary[[1L]]$lineage.coef[[lineage]]
  re <- primary[[2L]]$lineage.coef[[lineage]]
  if (!is.null(fe$error)) {
    stopifnot(!is.null(re$error), grepl("No studies available", fe$error))
    next
  }
  for (part in c("comp", "pa")) for (field in c("psi.lst", "cov.lst",
      "alignment", "coefficient.ids", "study.ids", "subject.ids")) {
    stopifnot(identical(fe[[paste0("score.summary.", part)]][[field]],
                        re[[paste0("score.summary.", part)]][[field]]))
  }
}

# Retain original local coordinates, dropping only globally unsupported columns.
component.bootstrap <- function(summary, keep.study, weights, method, component) {
  included <- summary$study.ids[summary$study.ids %in% keep.study]
  if (!length(included)) stop("No eligible studies remain for ", component, ".")
  psi <- summary$psi.lst[included]
  maps <- summary$alignment[included]
  supported <- sort(unique(unlist(maps, use.names = FALSE)))
  maps <- lapply(maps, function(x) match(x, supported))
  model <- list(cluster.score.lst = lapply(psi, t),
                bread.lst = lapply(psi, function(x) diag(ncol(x))),
                cov.lst = summary$cov.lst[included])
  multiplier <- lapply(included, function(study) {
    id <- rownames(psi[[study]])
    if (!all(id %in% rownames(weights[[study]]))) stop("Subject keys do not align.")
    weights[[study]][id, , drop = FALSE]
  })
  out <- MetaQCAT:::score_bootstrap_meta(model, multiplier, maps,
    length(supported), 1L, sub("-O$", "", method), component)
  draws <- ncol(multiplier[[1L]])
  stopifnot(length(out$bootstrap.statistic) == draws,
            all(is.finite(c(out$observed.statistic, out$bootstrap.statistic))))
  out$p <- (1 + sum(out$bootstrap.statistic >= out$observed.statistic)) / (draws + 1)
  out$study.beta.reps <- NULL
  out$included.studies <- included
  out$coefficient.ids <- summary$coefficient.ids[supported]
  out
}

# Recover every successful full-data statistic/draw/p-value before any omission.
verification.weights <- MetaQCAT:::.make_cluster_weights(subject_ids, 1000L, seed)
verification.rows <- list()
for (method in methods) for (lineage in lineages) {
  original <- primary[[method]]$lineage.coef[[lineage]]
  if (!is.null(original$error)) next
  marginal <- numeric(2)
  for (k in 1:2) {
    part <- c("comp", "pa")[k]
    recovered <- component.bootstrap(original[[paste0("score.summary.", part)]],
      names(counts), verification.weights, method, c("count", "zero")[k])
    marginal[k] <- recovered$p
    verification.rows[[paste(method, lineage, part)]] <- data.frame(
      method = method, lineage = lineage, component = c("count", "zero")[k],
      observed_error = abs(recovered$observed.statistic - original[[paste0("score.stat.", part)]]),
      max_bootstrap_error = max(abs(recovered$bootstrap.statistic - original[[paste0("bootstrap.stat.", part)]])),
      original_p = unname(original[[paste0("pvalue.", part)]]["Boot"]),
      recovered_p = recovered$p, stringsAsFactors = FALSE)
  }
  stopifnot(isTRUE(all.equal(MetaQCAT:::ACAT(marginal),
    unname(original$pvalue.comb["Boot"]), tolerance = 1e-14)))
}
verification <- do.call(rbind, verification.rows)
stopifnot(nrow(verification) == 168L, all(verification$observed_error == 0),
          all(verification$max_bootstrap_error == 0),
          identical(verification$original_p, verification$recovered_p))
write.csv(verification, file.path(result.dir,
  "IBD_primary_score_recovery_verification.csv"), row.names = FALSE)
rm(verification.weights)

# One shared 300-draw plan for every cohort omission, part, method, and lineage.
weights <- MetaQCAT:::.make_cluster_weights(subject_ids, n.boot, seed)
weights.sha256 <- digest::digest(weights, algo = "sha256")
lineage.rows <- run.rows <- list()
for (method in methods) {
  fit <- primary[[method]]
  p.row <- paste0("Bootstrap-", method)
  p <- fit$lineage.pval[p.row, lineages]
  q <- fit$lineage.qval[p.row, lineages]
  reason <- setNames(rep("", length(lineages)), lineages)
  reason[names(fit$failed.lineage)] <- fit$failed.lineage
  lineage.rows[[paste("Full_data_1000", method)]] <- data.frame(
    omitted_cohort = "Full_data_1000", method = method, lineage = lineages,
    p = unname(p), q = unname(q), significant = ifelse(is.na(q), NA, q < .05),
    failure_reason = unname(reason), stringsAsFactors = FALSE)
}

for (cohort in c("Full_data_300", cohorts)) {
  omitted.strata <- if (cohort == "Full_data_300") character() else
    names(source.cohort)[source.cohort == cohort]
  keep.study <- setdiff(names(counts), omitted.strata)
  message("IBD score LOCO: ", cohort, "; ", length(keep.study), " strata remain.")
  for (method in methods) {
    started <- Sys.time()
    p <- setNames(rep(NA_real_, length(lineages)), lineages)
    reason <- setNames(rep("", length(lineages)), lineages)
    detailed <- setNames(vector("list", length(lineages)), lineages)
    for (lineage in lineages) {
      original <- primary[[1L]]$lineage.coef[[lineage]]
      if (!is.null(original$error)) {
        reason[lineage] <- original$error
        detailed[[lineage]] <- list(error = original$error)
        next
      }
      result <- tryCatch({
        comp <- component.bootstrap(original$score.summary.comp, keep.study,
                                     weights, method, "count")
        zero <- component.bootstrap(original$score.summary.pa, keep.study,
                                     weights, method, "zero")
        list(count = comp, zero = zero, p = MetaQCAT:::ACAT(c(comp$p, zero$p)))
      }, error = function(e) e)
      if (inherits(result, "error")) {
        reason[lineage] <- conditionMessage(result)
        detailed[[lineage]] <- list(error = conditionMessage(result))
      } else {
        p[lineage] <- result$p
        detailed[[lineage]] <- result
      }
    }
    q <- p.adjust(p, method = "BH")
    fit <- list(lineage.pval = matrix(p, 1, dimnames = list(paste0("Bootstrap-", method), lineages)),
      lineage.qval = matrix(q, 1, dimnames = list(paste0("Bootstrap-", method), lineages)),
      lineage.coef = detailed, failed.lineage = reason[nzchar(reason)],
      bootstrap = list(calibration = "null-imposed wild-score bootstrap", requested = n.boot,
        package.version = as.character(utils::packageVersion("MetaQCAT")),
        re.variance.score.formula = "positive-part for count and zero",
        completed = if (any(nzchar(reason))) 0L else n.boot,
        completed.by.lineage = setNames(ifelse(is.na(p), 0L, n.boot), lineages),
        n.clusters = vapply(weights[keep.study], nrow, integer(1)),
        weights.sha256 = weights.sha256, null.models.refitted = FALSE,
        seed = seed, frozen.primary.preprocessing = TRUE,
        primary.file = primary.file, primary.md5 = source.md5))
    safe.cohort <- gsub("[^A-Za-z0-9_.-]", "_", cohort)
    saveRDS(fit, file.path(result.dir, paste0("score_", safe.cohort, "_",
      gsub("[^A-Za-z0-9]+", "_", method), ".rds")))
    lineage.rows[[paste(cohort, method)]] <- data.frame(omitted_cohort = cohort,
      method = method, lineage = lineages, p = unname(p), q = unname(q),
      significant = ifelse(is.na(q), NA, q < .05), failure_reason = unname(reason),
      stringsAsFactors = FALSE)
    run.rows[[paste(cohort, method)]] <- data.frame(omitted_cohort = cohort,
      method = method, omitted_strata = paste(omitted.strata, collapse = ";"),
      remaining_strata = length(keep.study),
      samples = sum(vapply(counts[keep.study], nrow, integer(1))),
      cases = sum(vapply(design[keep.study], function(x) sum(x[, "disease"] == 1), numeric(1))),
      controls = sum(vapply(design[keep.study], function(x) sum(x[, "disease"] == 0), numeric(1))),
      clusters = sum(vapply(subject_ids[keep.study], function(x) length(unique(x)), integer(1))),
      bootstrap_budget = n.boot, seed = seed, requested_threads = 1L,
      workers_used = 1L, seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      failed_lineages = sum(nzchar(reason)), error = "", frozen_primary_preprocessing = TRUE,
      stringsAsFactors = FALSE)
  }
}

all.lineages <- do.call(rbind, lineage.rows)
stability.results <- list()
for (baseline.label in c("Full_data_300", "Full_data_1000")) {
  stability.rows <- list()
  for (cohort in c(baseline.label, cohorts)) for (method in methods) {
    current <- all.lineages[all.lineages$omitted_cohort == cohort & all.lineages$method == method, ]
    baseline <- all.lineages[all.lineages$omitted_cohort == baseline.label & all.lineages$method == method, ]
    primary.sig <- baseline$lineage[!is.na(baseline$q) & baseline$q < .05]
    current.sig <- current$lineage[!is.na(current$q) & current$q < .05]
    retained <- intersect(primary.sig, current.sig)
    union.sig <- union(primary.sig, current.sig)
    stability.rows[[paste(cohort, method)]] <- data.frame(omitted_cohort = cohort,
      method = method, primary_discoveries = length(primary.sig), discoveries = length(current.sig),
      retained = length(retained), lost = length(setdiff(primary.sig, current.sig)),
      gained = length(setdiff(current.sig, primary.sig)),
      retention_fraction = length(retained) / length(primary.sig),
      jaccard = length(retained) / length(union.sig),
      retained_lineages = paste(retained, collapse = ";"),
      lost_lineages = paste(setdiff(primary.sig, current.sig), collapse = ";"),
      gained_lineages = paste(setdiff(current.sig, primary.sig), collapse = ";"),
      stringsAsFactors = FALSE)
  }
  stability.results[[baseline.label]] <- do.call(rbind, stability.rows)
}
budget.rows <- list()
for (method in methods) {
  a <- all.lineages[all.lineages$omitted_cohort == "Full_data_1000" & all.lineages$method == method, ]
  b <- all.lineages[all.lineages$omitted_cohort == "Full_data_300" & all.lineages$method == method, ]
  a <- a$lineage[!is.na(a$q) & a$q < .05]
  b <- b$lineage[!is.na(b$q) & b$q < .05]
  budget.rows[[method]] <- data.frame(method = method, discoveries_1000 = length(a),
    discoveries_300 = length(b), retained = length(intersect(a, b)),
    lost_at_300 = length(setdiff(a, b)), gained_at_300 = length(setdiff(b, a)),
    jaccard = length(intersect(a, b)) / length(union(a, b)),
    lost_lineages = paste(setdiff(a, b), collapse = ";"),
    gained_lineages = paste(setdiff(b, a), collapse = ";"))
}
write.csv(all.lineages, file.path(result.dir, "IBD_loso_lineage_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, run.rows), file.path(result.dir, "IBD_loso_run_summary.csv"), row.names = FALSE)
write.csv(stability.results$Full_data_300, file.path(result.dir, "IBD_loso_stability_summary.csv"), row.names = FALSE)
write.csv(stability.results$Full_data_1000, file.path(result.dir, "IBD_loco_primary_stability_summary.csv"), row.names = FALSE)
write.csv(stability.results$Full_data_300, file.path(result.dir, "IBD_loco_matched300_stability_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, budget.rows), file.path(result.dir, "IBD_loso_budget_comparison.csv"), row.names = FALSE)
saveRDS(list(settings = list(n.boot = n.boot, seed = seed, min.prev = .1,
  comparison_baseline = "Full_data_300",
  package.version = as.character(utils::packageVersion("MetaQCAT")),
  re.variance.score.formula = "positive-part for count and zero",
  weights.sha256 = weights.sha256, null.models.refitted = FALSE,
  source_primary = primary.file, source_primary_md5 = source.md5,
  source_data = data.file, source_data_md5 = unname(tools::md5sum(data.file)),
  frozen_primary_preprocessing = TRUE,
  preprocessing = "Primary rarefaction, per-study eligibility, and reference taxa retained; scores not refitted.",
  total_seconds = as.numeric(difftime(Sys.time(), script.started, units = "secs"))),
  lineage_results = all.lineages, run_summary = do.call(rbind, run.rows),
  stability_summary = stability.results$Full_data_300,
  primary_stability_summary = stability.results$Full_data_1000,
  matched300_stability_summary = stability.results$Full_data_300,
  budget_comparison = do.call(rbind, budget.rows), verification = verification,
  session_info = sessionInfo()), file.path(result.dir, "IBD_loso_results.rds"))
message("IBD frozen-score LOCO complete: ", result.dir)
