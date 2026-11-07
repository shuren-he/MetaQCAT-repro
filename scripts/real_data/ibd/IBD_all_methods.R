# IBD: shared paper preprocessing, all MetaQCAT tests, and QCAT-C.
# Run from the repository root, or execute the sections in RStudio.
# Running this file recomputes the analyses and replaces its result files.
source(file.path("R", "project_paths.R"), local = TRUE)

## 1. Paper settings ----------------------------------------------------------
min.depth <- 3000L             # retain raw sample depth strictly greater than 3000
max.zero.fraction <- 0.90
min.prev <- 0.10               # MetaQCAT component-specific prevalence filter
n.boot <- 1000L
n.perm <- 1000L                # QCAT-C between-subject label permutations
n.threads <- 7L
seed <- 1234L
fdr.alpha <- 0.05
for (package in c("MetaQCAT", "miLineage", "data.table", "dplyr")) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install the analysis dependencies first: missing ", package)
}
stopifnot(packageVersion("MetaQCAT") == "0.2.1")
data.dir <- file.path(project.dir, "data", "IBD")
result.dir <- file.path(project.dir, "res", "IBD")
bootstrap.dir <- file.path(result.dir, "wild_score_bootstrap")
table.dir <- file.path(project.dir, "tables")

## 2. Load the raw data ONCE ---------------------------------------------------
count.env <- new.env(parent = emptyenv())
meta.env <- new.env(parent = emptyenv())
tax.env <- new.env(parent = emptyenv())
load(file.path(data.dir, "count.genus.RData"), envir = count.env)
load(file.path(data.dir, "meta.RData"), envir = meta.env)
load(file.path(data.dir, "tax.RData"), envir = tax.env)
otu <- count.env$otu
raw.metadata <- as.data.frame(meta.env$meta.data, stringsAsFactors = FALSE)
taxonomy <- as.matrix(tax.env$Tax[, seq_len(6L), drop = FALSE])
colnames(taxonomy) <- paste0("Rank", seq_len(6L))
stopifnot(identical(colnames(otu), rownames(raw.metadata)),
          identical(rownames(otu), rownames(taxonomy)),
          !anyDuplicated(colnames(otu)), !anyNA(otu), all(otu >= 0))
depth <- colSums(otu)

# These two small helpers only record exclusions and select binary covariates.
record.exclusion <- function(m, reason) {
  data.frame(sample_id = rownames(m), source_cohort = m$dataset_name,
    depth = unname(depth[rownames(m)]), reason = rep(reason, nrow(m)),
    setting = rep("primary", nrow(m)),
    stringsAsFactors = FALSE)
}
binary.estimable <- function(x) {
  tab <- table(x[!is.na(x)])
  length(tab) == 2L && max(tab) / sum(tab) < 0.95
}

## 3. ONE shared sample-selection procedure -----------------------------------
m <- raw.metadata
exclusions <- list()
excluded.cohorts <- c("PROTECT", "LSS-PRISM", "BIDMC-FMT")
drop <- m$dataset_name %in% excluded.cohorts
exclusions[[length(exclusions) + 1L]] <- record.exclusion(
  m[drop, , drop = FALSE], "excluded_source_cohort")
m <- m[!drop, , drop = FALSE]
m$age.numeric <- suppressWarnings(as.numeric(as.character(m$age)))
drop <- m$dataset_name == "HMP2" & is.na(m$age.numeric)
exclusions[[length(exclusions) + 1L]] <- record.exclusion(
  m[drop, , drop = FALSE], "HMP2_age_missing")
m <- m[!drop, , drop = FALSE]
drop <- depth[rownames(m)] <= min.depth
exclusions[[length(exclusions) + 1L]] <- record.exclusion(
  m[drop, , drop = FALSE], "raw_depth_at_or_below_3000")
m <- m[!drop, , drop = FALSE]

# Analyze specimen types separately and split HMP2 at age 18.
m$analysis_study <- m$dataset_name
hmp2 <- m$dataset_name == "HMP2"
m$analysis_study[hmp2 & m$age.numeric < 18] <- "HMP2.pediatric"
m$analysis_study[hmp2 & m$age.numeric >= 18] <- "HMP2.adult"
m$study_key <- paste(m$analysis_study, m$sample_type, sep = ".")
stopifnot(!anyNA(m$study_key))

# A subject cannot occur in more than one analysis stratum. Keep biopsy
# observations when that subject also contributes stool observations.
subject.key <- paste(m$dataset_name, m$subject_accession, sep = "::")
cross.stratum <- vapply(split(m$study_key, subject.key),
  function(x) length(unique(x)) > 1L, logical(1))
drop <- subject.key %in% names(cross.stratum)[cross.stratum] &
  m$sample_type == "stool"
exclusions[[length(exclusions) + 1L]] <- record.exclusion(
  m[drop, , drop = FALSE], "cross_stratum_overlap_keep_biopsy")
m <- m[!drop, , drop = FALSE]
subject.key <- paste(m$dataset_name, m$subject_accession, sep = "::")
stopifnot(!any(vapply(split(m$study_key, subject.key),
  function(x) length(unique(x)) > 1L, logical(1))))

## 4. Shared genus screening and adjusted design matrices ----------------------
by.study <- split(m, m$study_key, drop = TRUE)
counts <- metadata <- design <- subject_ids <- list()
study.rows <- list()
for (study in names(by.study)) {
  sm <- by.study[[study]]
  sex <- tolower(trimws(as.character(sm$gender)))
  sex[!sex %in% c("f", "m")] <- NA_character_
  race <- tolower(trimws(as.character(sm$race)))
  race[is.na(race) | !nzchar(race)] <- NA_character_
  sex.male <- ifelse(is.na(sex), NA_real_, as.numeric(sex == "m"))
  race.nonwhite <- ifelse(is.na(race), NA_real_, as.numeric(race != "white"))
  age.years <- sm$age.numeric
  covariates <- list()
  if (mean(is.na(sex.male)) < 0.5 && binary.estimable(sex.male))
    covariates$sex_male <- sex.male
  if (mean(is.na(race.nonwhite)) < 0.5 && binary.estimable(race.nonwhite))
    covariates$race_nonwhite <- race.nonwhite
  if (mean(is.na(age.years)) < 0.5 &&
      length(unique(age.years[!is.na(age.years)])) > 1L)
    covariates$age_10y <- age.years

  # Complete cases only for the selected adjustment covariates. No carry-forward
  # imputation and no missing-indicator alternative are used in this analysis.
  if (length(covariates)) {
    keep <- Reduce(`&`, lapply(covariates, function(x) !is.na(x)))
    if (any(!keep)) exclusions[[length(exclusions) + 1L]] <- record.exclusion(
      sm[!keep, , drop = FALSE], "selected_adjustment_covariate_missing")
    sm <- sm[keep, , drop = FALSE]
    covariates <- lapply(covariates, function(x) x[keep])
  }
  if (!nrow(sm)) next
  disease.status <- as.integer(sm$disease != "control")
  if (anyNA(disease.status) || length(unique(disease.status)) < 2L)
    stop("Both disease groups are required in ", study)
  X <- matrix(disease.status, ncol = 1L,
    dimnames = list(rownames(sm), "disease"))
  for (column in names(covariates)) {
    value <- covariates[[column]]
    if (column == "age_10y") value <- (value - mean(value)) / 10
    candidate <- cbind(X, value)
    if (qr(cbind(1, candidate))$rank == ncol(candidate) + 1L) {
      X <- candidate
      colnames(X)[ncol(X)] <- column
    }
  }
  storage.mode(X) <- "double"
  y <- t(otu[, rownames(sm), drop = FALSE])
  y <- y[, colMeans(y == 0) <= max.zero.fraction, drop = FALSE]
  if (!ncol(y)) stop("No genera retained in ", study)
  stopifnot(identical(rownames(y), rownames(X)), !anyNA(X),
            qr(cbind(1, X))$rank == ncol(X) + 1L)
  counts[[study]] <- y
  metadata[[study]] <- sm
  design[[study]] <- X
  subject_ids[[study]] <- factor(sm$subject_accession)
  study.rows[[study]] <- data.frame(
    setting = "primary", study_key = study, source_cohort = sm$dataset_name[1L],
    specimen = sm$sample_type[1L], samples = nrow(sm),
    subjects = nlevels(subject_ids[[study]]),
    cases = sum(disease.status == 1L), controls = sum(disease.status == 0L),
    low_control_flag = sum(disease.status == 0L) < 10L,
    taxa = ncol(y), min_depth = min(rowSums(y)),
    median_depth = median(rowSums(y)),
    covariates = paste(colnames(X)[-1L], collapse = ";"),
    stringsAsFactors = FALSE)
}
study.summary <- do.call(rbind, study.rows)
rownames(study.summary) <- NULL
exclusions <- do.call(rbind, exclusions)
retained.ids <- unlist(lapply(metadata, rownames), use.names = FALSE)
stopifnot(length(counts) >= 2L, !anyDuplicated(retained.ids),
  !any(retained.ids %in% exclusions$sample_id),
  length(retained.ids) + length(unique(exclusions$sample_id)) == nrow(raw.metadata),
  identical(names(counts), names(design)),
  identical(names(counts), names(subject_ids)))
print(study.summary)

## 5. Save the ONE prepared dataset used by EVERY method ------------------------
dir.create(bootstrap.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table.dir, recursive = TRUE, showWarnings = FALSE)
disease <- lapply(design, function(x) x[, "disease", drop = FALSE])
adjust <- lapply(design, function(x) x[, -1L, drop = FALSE])
save(counts, metadata, design, subject_ids, taxonomy, disease, adjust,
  file = file.path(bootstrap.dir, "IBD_processed_bootstrap.RData"), compress = "xz")
write.csv(exclusions, file.path(bootstrap.dir, "IBD_nonoverlap_exclusions.csv"),
          row.names = FALSE)
write.csv(study.summary, file.path(bootstrap.dir, "IBD_nonoverlap_study_summary.csv"),
          row.names = FALSE)
write.csv(study.summary, file.path(bootstrap.dir, "IBD_setting_study_summary.csv"),
          row.names = FALSE)
settings <- list(phenotype = "IBD", specimen = "separate strata",
  min_raw_depth = min.depth, max_zero_fraction = max.zero.fraction,
  min.prev = min.prev, overlap_priority = "biopsy", missing_policy = "complete_case",
  adjusted = TRUE, n.boot = n.boot, n.perm = n.perm,
  n.threads = n.threads, seed = seed,
  input_md5 = tools::md5sum(file.path(data.dir,
    c("count.genus.RData", "meta.RData", "tax.RData"))),
  script_md5 = unname(tools::md5sum(file.path(project.dir, "scripts", "real_data",
                                             "ibd", "IBD_all_methods.R"))))

# Atomic writes protect each completed fit if a later analysis is interrupted.
save.result <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path)), tmpdir = dirname(path))
  saveRDS(x, tmp, compress = "gzip")
  if (!file.rename(tmp, path)) stop("Could not save result: ", path)
}

## 6. MetaQCAT: two fits provide all six component/omnibus tests ------------------
tax <- taxonomy[, c("Rank5", "Rank6"), drop = FALSE]
colnames(tax) <- c("Rank1", "Rank2")
run.started <- Sys.time()
started <- Sys.time()
FE_IBD <- MetaQCAT::MetaQCAT(
  counts = counts, X = design, X.index = 1L, cluster.id = subject_ids,
  taxonomy = tax, method = "FE-MetaQCAT-O", min.depth = 0, min.prev = min.prev,
  n.boot = n.boot, n.threads = n.threads, seed = seed)
fe.seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
save.result(FE_IBD, file.path(bootstrap.dir, "IBD_primary_FE_MetaQCAT_O.rds"))
started <- Sys.time()
RE_IBD <- MetaQCAT::MetaQCAT(
  counts = counts, X = design, X.index = 1L, cluster.id = subject_ids,
  taxonomy = tax, method = "RE-MetaQCAT-O", min.depth = 0, min.prev = min.prev,
  n.boot = n.boot, n.threads = n.threads, seed = seed)
re.seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
save.result(RE_IBD, file.path(bootstrap.dir, "IBD_primary_RE_MetaQCAT_O.rds"))
results <- setNames(list(FE_IBD, RE_IBD), c("FE-MetaQCAT-O", "RE-MetaQCAT-O"))
save.result(list(setting = settings, data_summary = study.summary, results = results),
            file.path(bootstrap.dir, "IBD_primary.rds"))
save.result(list(primary = results), file.path(bootstrap.dir, "IBD_all_results.rds"))

## 7. QCAT-C: the SAME counts, full design matrices, and subject clusters ---------
library(data.table)
source(file.path(project.dir, "code", "IBD", "QCAT_C.R"), local = TRUE)
old.rng.kind <- RNGkind()
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
started <- Sys.time()
QCAT_C_IBD <- tryCatch(QCAT_C_Meta(
  ID = subject_ids, OTU = counts, X = design, X.index = 1L, Tax = tax,
  perm.type = "BTW", n.perm = n.perm, test = "chisq", nthreads = n.threads
), error = identity, finally = do.call(RNGkind, as.list(old.rng.kind)))
qc.seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
save.result(QCAT_C_IBD, file.path(result.dir, "IBD_primary_QCAT_C.rds"))
if (inherits(QCAT_C_IBD, "error")) stop(conditionMessage(QCAT_C_IBD))
stopifnot("Resampling-Fisher" %in% rownames(QCAT_C_IBD$lineage.pval))
qc.p <- QCAT_C_IBD$lineage.pval["Resampling-Fisher", ]
qc.q <- p.adjust(qc.p, method = "BH")

## 8. Fresh p/q tables: all methods, no old result templates ----------------------
lineages <- colnames(FE_IBD$lineage.pval)
stopifnot(setequal(lineages, colnames(RE_IBD$lineage.pval)),
          !anyDuplicated(lineages))
method.names <- c(rep(c("FE-MetaQCAT", "FE-MetaQCAT*", "FE-MetaQCAT-O"), 2L),
                  "RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O", "QCAT-C")
calibrations <- c(rep("Asy", 3L), rep("Boot", 6L), "Permutation")
fields <- rep(c("pvalue.comp", "pvalue.pa", "pvalue.comb"), 3L)
p.values <- matrix(NA_real_, length(lineages), 10L,
  dimnames = list(lineages, c(paste0(method.names[1:9], " (",
    c(rep("asy", 3L), rep("boot", 6L)), ")"), "QCAT-C")))
for (j in seq_len(9L)) {
  fit <- if (j <= 6L) FE_IBD else RE_IBD
  p.values[, j] <- vapply(fit$lineage.coef[lineages], function(z) {
    p <- z[[fields[j]]]
    if (is.null(p) || !calibrations[j] %in% names(p)) NA_real_ else
      unname(p[calibrations[j]])
  }, numeric(1))
}
q.values <- apply(p.values, 2L, p.adjust, method = "BH")
# QCAT-C q-values are calculated over its own tested lineages first, then
# aligned to the MetaQCAT lineage universe; absent lineages remain NA.
p.values[, 10L] <- qc.p[match(lineages, names(qc.p))]
q.values[, 10L] <- qc.q[match(lineages, names(qc.q))]
stopifnot(isTRUE(all.equal(unname(q.values[, 6L]),
  unname(FE_IBD$lineage.qval["Bootstrap-FE-MetaQCAT-O", lineages]))),
  isTRUE(all.equal(unname(q.values[, 9L]),
  unname(RE_IBD$lineage.qval["Bootstrap-RE-MetaQCAT-O", lineages]))))
write.csv(p.values, file.path(result.dir, "IBD_family_p_val.csv"))
write.csv(q.values, file.path(table.dir, "IBD_family_q_val.csv"))
write.csv(q.values[, 4:10, drop = FALSE],
          file.path(table.dir, "IBD_family_q_val_resampling.csv"))

lineage.rows <- list()
for (j in c(6L, 9L)) {
  component.cols <- if (j == 6L) 4:5 else 7:8
  lineage.rows[[method.names[j]]] <- data.frame(
    setting = "primary", method = method.names[j], lineage = lineages,
    p = p.values[, j], q = q.values[, j],
    positive_p = p.values[, component.cols[1L]],
    positive_q = q.values[, component.cols[1L]],
    zero_p = p.values[, component.cols[2L]],
    zero_q = q.values[, component.cols[2L]],
    zero_only = !is.na(q.values[, j]) & q.values[, j] < fdr.alpha &
      !is.na(q.values[, component.cols[1L]]) &
      q.values[, component.cols[1L]] >= fdr.alpha &
      !is.na(q.values[, component.cols[2L]]) &
      q.values[, component.cols[2L]] < fdr.alpha,
    row.names = NULL, stringsAsFactors = FALSE)
}
latest <- do.call(rbind, lineage.rows)
rownames(latest) <- NULL
write.csv(latest, file.path(bootstrap.dir, "IBD_lineage_results.csv"), row.names = FALSE)
write.csv(latest, file.path(result.dir, "IBD_lineage_results.csv"), row.names = FALSE)
qc.table <- data.frame(setting = "primary", method = "QCAT-C", lineage = lineages,
  p = p.values[, 10L], q = q.values[, 10L], positive_p = NA_real_,
  positive_q = NA_real_, zero_p = NA_real_, zero_q = NA_real_,
  zero_only = FALSE, row.names = NULL, stringsAsFactors = FALSE)
write.csv(qc.table, file.path(result.dir, "IBD_QCAT_C_lineage_results.csv"),
          row.names = FALSE)

# Match the existing UpSet/extraction interface without importing old tables.
fe <- lineage.rows[["FE-MetaQCAT-O"]]
re <- lineage.rows[["RE-MetaQCAT-O"]]
combined <- data.frame(lineage = lineages, gap1 = NA_character_,
  method = "FE-MetaQCAT-O", p = fe$p, q = fe$q, zero_only = fe$zero_only,
  gap2 = NA_character_, method.1 = "RE-MetaQCAT-O", p.1 = re$p, q.1 = re$q,
  zero_only.1 = re$zero_only, gap3 = NA_character_, method.2 = "QCAT-C",
  p.2 = qc.table$p, q.2 = qc.table$q, check.names = FALSE)
names(combined)[c(2L, 7L, 12L)] <- ""
write.csv(combined, file.path(result.dir, "IBD_lineage_combine.csv"), row.names = FALSE)
sig <- data.frame(Lineage = lineages, check.names = FALSE)
sig[["QCAT-C"]] <- as.integer(!is.na(qc.table$q) & qc.table$q < fdr.alpha)
sig$FE <- as.integer(!is.na(fe$q) & fe$q < fdr.alpha)
sig$RE <- as.integer(!is.na(re$q) & re$q < fdr.alpha)
write.csv(sig, file.path(result.dir, "IBD_sig_family.csv"), row.names = FALSE)
discovery.summary <- data.frame(
  strategy = c("QCAT-C", "FE", "RE", "FE_and_RE", "FE_only", "RE_only", "Any"),
  n_discoveries = c(colSums(sig[, -1L]), sum(sig$FE == 1L & sig$RE == 1L),
    sum(sig$FE == 1L & sig$RE == 0L), sum(sig$FE == 0L & sig$RE == 1L),
    sum(rowSums(sig[, -1L]) > 0L)))
write.csv(discovery.summary, file.path(bootstrap.dir, "discovery_summary.csv"),
          row.names = FALSE)

## 9. Save complete method results and timing -----------------------------------
all.method.rows <- lapply(seq_along(method.names), function(j) data.frame(
  lineage = lineages, method = method.names[j], calibration = calibrations[j],
  p = p.values[, j], q = q.values[, j], row.names = NULL))
write.csv(do.call(rbind, all.method.rows),
          file.path(result.dir, "IBD_all_method_results.csv"), row.names = FALSE)
run.summary <- data.frame(setting = "primary",
  method = c("FE-MetaQCAT-O", "RE-MetaQCAT-O", "QCAT-C"),
  samples = sum(vapply(counts, nrow, integer(1))), strata = length(counts),
  source_cohorts = length(unique(study.summary$source_cohort)),
  bootstrap_budget = c(n.boot, n.boot, 0L), permutation_budget = c(0L, 0L, n.perm),
  seed = seed, seconds = c(fe.seconds, re.seconds, qc.seconds),
  failed_lineages = c(length(FE_IBD$failed.lineage), length(RE_IBD$failed.lineage),
                       sum(is.na(qc.table$p))),
  n.threads = n.threads, stringsAsFactors = FALSE)
write.csv(run.summary, file.path(bootstrap.dir, "IBD_run_summary.csv"), row.names = FALSE)
save.result(list(settings = settings, results = results, QCAT_C = QCAT_C_IBD,
  p_values = p.values, q_values = q.values, significant = sig,
  study_summary = study.summary, exclusions = exclusions, run_summary = run.summary,
  elapsed_seconds = as.numeric(difftime(Sys.time(), run.started, units = "secs")),
  session_info = sessionInfo()), file.path(result.dir, "IBD_all_methods.rds"))
message("Saved IBD MetaQCAT and QCAT-C analyses under ", result.dir)
