# IBD: current paper preprocessing, MetaQCAT, and cluster QCAT under P/CD/RD.
# Run from the repository root, or execute the sections in RStudio.
# Running this file recomputes the analyses and replaces its result files.
source(file.path("R", "project_paths.R"), local = TRUE)

## 1. Paper settings ----------------------------------------------------------
min.depth <- 3000L             # retain raw sample depth strictly greater than 3000
max.zero.fraction <- 0.90
min.prev <- 0.10               # MetaQCAT component-specific prevalence filter
n.boot <- 1000L
n.perm <- as.integer(Sys.getenv("IBD_QCAT_PERMS", "300")) # native BTW cap
n.threads <- as.integer(Sys.getenv("IBD_QCAT_CORES", "7"))
seed <- 1234L
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
fdr.alpha <- 0.05
run.metaqcat <- match.arg(tolower(Sys.getenv("IBD_RUN_METAQCAT", "true")),
                         c("true", "false")) == "true"
prepare.only <- match.arg(tolower(Sys.getenv("IBD_PREPARE_ONLY", "false")),
                         c("true", "false")) == "true"
stopifnot(!anyNA(c(n.perm, n.threads)), n.perm >= 1L, n.threads >= 1L)
if (.Platform$OS.type == "windows") n.threads <- 1L
available.cores <- parallel::detectCores()
if (is.finite(available.cores)) n.threads <- min(n.threads, available.cores)
for (package in c("MetaQCAT", "miLineage", "MMUPHin", "digest")) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install the analysis dependencies first: missing ", package)
}
stopifnot(packageVersion("MetaQCAT") == "0.2.1")
data.dir <- file.path(project.dir, "data", "IBD")
if (!dir.exists(data.dir))
  data.dir <- file.path(project.dir, "real", "IBD", "data")
if (!all(file.exists(file.path(data.dir,
  c("count.genus.RData", "meta.RData", "tax.RData")))))
  stop("Place the three IBD inputs under data/IBD/ or real/IBD/data/.")
result.dir <- file.path(project.dir, "res", "IBD")
bootstrap.dir <- file.path(result.dir, "wild_score_bootstrap")
table.dir <- file.path(project.dir, "tables")
cluster.dir <- file.path(result.dir, "cluster_qcat")
prepared.file <- file.path(bootstrap.dir, "IBD_processed_bootstrap.RData")
all.methods.file <- file.path(project.dir, "scripts", "real_data", "ibd",
                              "IBD_all_methods.R")

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

## 5. Save and audit the shared preprocessing ---------------------------------
# Reuse of old FE/RE fits is allowed only when their saved analysis inputs agree.
if (!run.metaqcat && !prepare.only) {
  if (!file.exists(prepared.file)) stop("No saved IBD primary preprocessing to reuse.")
  prior <- new.env(parent = emptyenv())
  load(prepared.file, envir = prior)
  for (object in c("counts", "metadata", "design", "subject_ids", "taxonomy"))
    if (!identical(get(object), get(object, envir = prior, inherits = FALSE)))
      stop("Saved primary inputs differ: ", object, ". Rerun MetaQCAT first.")
}
dir.create(bootstrap.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table.dir, recursive = TRUE, showWarnings = FALSE)
disease <- lapply(design, function(x) x[, "disease", drop = FALSE])
adjust <- lapply(design, function(x) x[, -1L, drop = FALSE])
save(counts, metadata, design, subject_ids, taxonomy, disease, adjust,
  file = prepared.file, compress = "xz")
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
  perm.type = "BTW", run.metaqcat = run.metaqcat,
  n.threads = n.threads, seed = seed, rng = RNGkind(),
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

## 6. MetaQCAT: two fits provide all six component/omnibus tests ----------------
saveRDS(list(settings = settings, session_info = sessionInfo()),
        file.path(bootstrap.dir, "IBD_workflow_configuration.rds"))
if (prepare.only) {
  message("Preprocessing only: ", sum(vapply(counts, nrow, integer(1))),
          " samples; ", length(counts), " strata; ",
          sum(vapply(subject_ids, nlevels, integer(1))), " subject clusters.")
  if (interactive()) stop("Preprocessing complete; no models run.", call. = FALSE)
  quit(save = "no", status = 0L)
}
tax <- taxonomy[, c("Rank5", "Rank6"), drop = FALSE]
colnames(tax) <- c("Rank1", "Rank2")
run.started <- Sys.time()
fe.seconds <- re.seconds <- 0
if (run.metaqcat) {
  started <- Sys.time()
  FE_IBD <- MetaQCAT::MetaQCAT(
    counts = counts, X = design, X.index = 1L, cluster.id = subject_ids,
    taxonomy = tax, method = "FE-MetaQCAT-O", min.depth = 0, min.prev = min.prev,
    n.boot = n.boot, n.threads = n.threads, seed = seed)
  fe.seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  FE_IBD$bootstrap$package.version <- as.character(packageVersion("MetaQCAT"))
  FE_IBD$bootstrap$re.variance.score.formula <- "positive-part for count and zero"
  save.result(FE_IBD, file.path(bootstrap.dir, "IBD_primary_FE_MetaQCAT_O.rds"))
  started <- Sys.time()
  RE_IBD <- MetaQCAT::MetaQCAT(
    counts = counts, X = design, X.index = 1L, cluster.id = subject_ids,
    taxonomy = tax, method = "RE-MetaQCAT-O", min.depth = 0, min.prev = min.prev,
    n.boot = n.boot, n.threads = n.threads, seed = seed)
  re.seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  RE_IBD$bootstrap$package.version <- as.character(packageVersion("MetaQCAT"))
  RE_IBD$bootstrap$re.variance.score.formula <- "positive-part for count and zero"
  save.result(RE_IBD, file.path(bootstrap.dir, "IBD_primary_RE_MetaQCAT_O.rds"))
} else {
  FE_IBD <- readRDS(file.path(bootstrap.dir, "IBD_primary_FE_MetaQCAT_O.rds"))
  RE_IBD <- readRDS(file.path(bootstrap.dir, "IBD_primary_RE_MetaQCAT_O.rds"))
  for (fit in list(FE_IBD, RE_IBD)) {
    stopifnot(is.list(fit),
      identical(as.integer(fit$bootstrap$requested), as.integer(n.boot)),
      identical(as.integer(fit$bootstrap$seed), as.integer(seed)),
      identical(fit$bootstrap$package.version,
                as.character(packageVersion("MetaQCAT"))),
      setequal(names(fit$bootstrap$n.clusters), names(subject_ids)),
      identical(as.integer(fit$bootstrap$n.clusters[names(subject_ids)]),
                unname(vapply(subject_ids, nlevels, integer(1)))),
      identical(fit$bootstrap$re.variance.score.formula,
                "positive-part for count and zero"))
  }
  message("Reusing validated FE/RE primary fits; no MetaQCAT models rerun.")
}
meta.results <- setNames(list(FE_IBD, RE_IBD),
                         c("FE-MetaQCAT-O", "RE-MetaQCAT-O"))
lineages <- colnames(FE_IBD$lineage.pval)
stopifnot(setequal(lineages, colnames(RE_IBD$lineage.pval)),
          !anyDuplicated(lineages))
save.result(list(setting = settings, data_summary = study.summary,
                 results = meta.results), file.path(bootstrap.dir, "IBD_primary.rds"))
save.result(list(primary = meta.results),
            file.path(bootstrap.dir, "IBD_all_results.rds"))

## 7. Cluster QCAT comparisons: shared samples, native BTW permutation ----------
# QCAT1 below calls QCAT.Cluster (QCAT in the paper).
# QCAT2 calls QCAT_GEE.Cluster; its two-part input is rarefied once.
comparison.packages <- c("miLineage", "MMUPHin", "MetaQCAT", "digest")
families <- lineages
stopifnot("perm.type" %in% names(formals(miLineage::QCAT.Cluster)),
          "perm.type" %in% names(formals(miLineage::QCAT_GEE.Cluster)))
for (study in names(counts)) {
  status <- split(design[[study]][, "disease"], subject_ids[[study]])
  if (!all(vapply(status, function(z) length(unique(z)) == 1L, logical(1))))
    stop("Between-subject permutation requires subject-constant disease: ", study)
}
dir.create(file.path(cluster.dir, "method_results"),
           recursive = TRUE, showWarnings = FALSE)

# Pooled original counts on exactly the same retained samples.
# As in pour_methods.R, pooled strategies screen the original pooled table,
# not zero-filled versions of the separately screened study tables.
sample.ids <- unlist(lapply(counts, rownames), use.names = FALSE)
pooled <- t(otu[, sample.ids, drop = FALSE])
pooled <- pooled[, colMeans(pooled == 0) <= max.zero.fraction, drop = FALSE]
pooled.study <- rep(names(counts), vapply(counts, nrow, integer(1)))
pooled.id <- unlist(lapply(names(counts), function(study)
  paste(study, as.character(subject_ids[[study]]), sep = "::")), use.names = FALSE)
pooled.design <- matrix(unlist(lapply(design, function(x) x[, "disease"]),
                              use.names = FALSE), ncol = 1L,
  dimnames = list(sample.ids, "disease"))
# Selected adjustment sets differ by stratum. Each selected adjustment column
# is retained only for its own stratum (zero outside it), without imputing an
# unselected/missing covariate for another stratum. No study indicators are added.
for (study in names(design)) {
  rows <- which(pooled.study == study)
  for (column in setdiff(colnames(design[[study]]), "disease")) {
    value <- numeric(nrow(pooled.design))
    value[rows] <- design[[study]][, column]
    candidate <- cbind(pooled.design, value)
    if (qr(cbind(1, candidate))$rank == ncol(candidate) + 1L) {
      pooled.design <- candidate
      colnames(pooled.design)[ncol(pooled.design)] <-
        make.names(paste(study, column, sep = "__"))
    }
  }
}
stopifnot(identical(rownames(pooled), rownames(pooled.design)),
          all(rowSums(pooled) > 0), !anyNA(pooled.design))

# One MMUPHin correction and rarefaction before family tests.
# Checkpointed preprocessing is shared by all families and QCAT variants.
configuration <- list(n.perm = n.perm, seed = seed,
  rng = RNGkind(),
  rarefaction_seed = seed, native_permutation_seed = 11L,
  min_raw_depth = min.depth, max_zero_fraction = max.zero.fraction,
  perm.type = "BTW", test = "chisq", families = families,
  pooled_adjustment = "study-specific selected covariates; no study indicators",
  input.md5 = tools::md5sum(c(prepared.file,
    file.path(data.dir, "count.genus.RData"))),
  source.md5 = unname(tools::md5sum(all.methods.file)),
  versions = vapply(comparison.packages, function(p) as.character(utils::packageVersion(p)), character(1)))
fingerprint <- digest::digest(configuration)
saveRDS(list(configuration = configuration, n.threads = n.threads,
             session = sessionInfo()), file.path(cluster.dir, "configuration.rds"))
preprocessing.file <- file.path(cluster.dir, "preprocessing.rds")
cached <- if (file.exists(preprocessing.file)) readRDS(preprocessing.file) else NULL
if (!is.null(cached) && identical(cached$fingerprint, fingerprint)) {
  corrected <- cached$corrected
  rare.study <- cached$rare.study
  rare.pooled <- cached$rare.pooled
  rare.corrected <- cached$rare.corrected
  preprocessing.seconds <- cached$seconds
} else {
  started <- proc.time()[["elapsed"]]
  batch.metadata <- as.data.frame(pooled.design)
  batch.metadata$study <- factor(pooled.study)
  rownames(batch.metadata) <- sample.ids
  set.seed(seed)
  batch.fit <- MMUPHin::adjust_batch(
    feature_abd = t(pooled), batch = "study",
    covariates = colnames(pooled.design), data = batch.metadata,
    control = list(verbose = FALSE, diagnostic_plot = NULL))
  corrected <- t(batch.fit$feature_abd_adj)
  stopifnot(identical(dimnames(corrected), dimnames(pooled)),
    !anyNA(corrected), all(corrected >= 0),
    all(abs(corrected - round(corrected)) < 1e-8))
  corrected <- round(corrected)
  # Sequential hypergeometric rarefaction has the same sampling distribution
  # as the read-expansion Rarefy in the original comparison scripts.
  rarefy <- getFromNamespace(".rarefy_counts", "MetaQCAT")
  set.seed(seed)
  rare.study <- lapply(counts, rarefy)
  set.seed(seed)
  rare.pooled <- rarefy(pooled)
  set.seed(seed)
  rare.corrected <- rarefy(corrected)
  preprocessing.seconds <- proc.time()[["elapsed"]] - started
  saveRDS(list(fingerprint = fingerprint, corrected = corrected,
    rare.study = rare.study, rare.pooled = rare.pooled,
    rare.corrected = rare.corrected, seconds = preprocessing.seconds), preprocessing.file)
}

# Native cluster calls; miLineage functions are not modified.
run.cluster.test <- function(Y, X, ID, qcat2 = FALSE) {
  started <- proc.time()[["elapsed"]]
  warning.text <- character()
  error.text <- ""
  p <- c(Asymptotic = NA_real_, Resampling = NA_real_)
  x.names <- z.names <- character()
  count.disabled <- FALSE
  samples <- nrow(Y)
  clusters <- length(unique(ID))
  printed <- capture.output(answer <- tryCatch(withCallingHandlers({
    Y <- as.matrix(Y)
    Y <- Y[, colSums(Y) > 0, drop = FALSE]
    if (ncol(Y) < 2L) stop("Fewer than two observed genera.")
    stopifnot(nrow(Y) == nrow(X), length(ID) == nrow(Y),
      !anyNA(ID), !anyNA(X), all(rowSums(Y) >= 0))
    # Native BTW concatenates clusters; all its inputs must share this order.
    order.id <- order(as.character(ID), rownames(Y), method = "radix")
    Y <- Y[order.id, , drop = FALSE]
    X <- X[order.id, , drop = FALSE]
    ID <- factor(ID[order.id])
    status <- split(X[, 1L], ID)
    if (!all(vapply(status, function(z) length(unique(z)) == 1L, logical(1))))
      stop("Tested disease covariate varies within a subject.")
    if (length(unique(X[, 1L])) < 2L) stop("Only one disease group remains.")
    positive <- rowSums(Y) > 0
    # Select only an estimable nuisance basis. This leaves the nuisance span
    # unchanged and avoids native GEE's double-indexing of constant columns.
    nuisance.x <- nuisance.z <- integer()
    for (column in setdiff(seq_len(ncol(X)), 1L)) {
      x.cols <- c(nuisance.x, column)
      if (qr(cbind(1, X[positive, x.cols, drop = FALSE]))$rank == length(x.cols) + 1L)
        nuisance.x <- x.cols
      z.cols <- c(nuisance.z, column)
      if (qr(cbind(1, X[, z.cols, drop = FALSE]))$rank == length(z.cols) + 1L)
        nuisance.z <- z.cols
    }
    count.estimable <- sum(positive) > length(nuisance.x) + 2L &&
      qr(cbind(1, X[positive, c(1L, nuisance.x), drop = FALSE]))$rank ==
        length(nuisance.x) + 2L
    zero.estimable <- qr(cbind(1, X[, c(1L, nuisance.z), drop = FALSE]))$rank ==
      length(nuisance.z) + 2L
    if (!qcat2 && !count.estimable) stop("Disease is not estimable in positive-depth rows.")
    if (qcat2 && !zero.estimable) stop("Disease is aliased with the zero-part nuisance design.")
    XX <- X[, c(1L, nuisance.x), drop = FALSE]
    ZZ <- X[, c(1L, nuisance.z), drop = FALSE]
    x.names <- colnames(XX)
    z.names <- colnames(ZZ)
    count.disabled <- qcat2 && !count.estimable
    if (qcat2) {
      miLineage::QCAT_GEE.Cluster(ID = ID, OTU = Y, X = XX,
        X.index = if (count.disabled) integer(0) else 1L,
        Z = ZZ, Z.index = 1L, Tax = NULL, min.depth = 0,
        perm.type = "BTW", n.perm = n.perm, test = "chisq")
    } else {
      miLineage::QCAT.Cluster(ID = ID, OTU = Y, X = XX, X.index = 1L,
        Tax = NULL, min.depth = 0, perm.type = "BTW", n.perm = n.perm, test = "chisq")
    }
  }, warning = function(w) {
    warning.text <<- c(warning.text, conditionMessage(w))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")), error = identity))
  if (inherits(answer, "error")) error.text <- conditionMessage(answer) else {
    value <- as.numeric(answer$pval)
    if (length(value) >= 1L) p[1L] <- value[1L]
    if (length(value) >= 2L) p[2L] <- value[2L]
    p[!is.finite(p) | p < 0 | p > 1] <- NA_real_
    if (anyNA(p)) error.text <- "Native test returned a nonfinite p-value."
  }
  if (nzchar(error.text) && length(printed))
    error.text <- paste(error.text, substr(paste(tail(printed, 3L), collapse = " | "), 1L, 1200L), sep = " | ")
  list(p = p, error = error.text, warning = paste(unique(warning.text), collapse = " | "),
    seconds = proc.time()[["elapsed"]] - started, samples = samples, clusters = clusters,
    count_covariates = paste(x.names, collapse = ";"),
    zero_covariates = if (qcat2) paste(z.names, collapse = ";") else "",
    count_disabled = count.disabled)
}

# P combines studies; CD/RD analyze pooled data. Family checkpoints are resumable.
run.family <- function(family) {
  cache.file <- file.path(cluster.dir, "method_results",
    paste0(gsub("[^[:alnum:]_.-]", "_", family), ".rds"))
  previous <- if (file.exists(cache.file)) readRDS(cache.file) else NULL
  if (!is.null(previous) && identical(previous$fingerprint, fingerprint)) return(previous)
  started <- proc.time()[["elapsed"]]
  genera <- rownames(taxonomy)[taxonomy[, "Rank5"] == family]
  records <- list()
  calls <- list()
  for (strategy in c("P", "CD", "RD")) for (test in c("QCAT1", "QCAT2")) {
    qcat2 <- test == "QCAT2"
    if (strategy == "P") {
      studies <- names(counts)
      values <- matrix(NA_real_, length(studies), 2L,
        dimnames = list(studies, c("Asymptotic", "Resampling")))
      for (study in studies) {
        Y <- if (qcat2) rare.study[[study]] else counts[[study]]
        Y <- Y[, intersect(colnames(Y), genera), drop = FALSE]
        fit <- run.cluster.test(Y, design[[study]], subject_ids[[study]], qcat2)
        values[study, ] <- fit$p
        calls[[length(calls) + 1L]] <- data.frame(lineage = family,
          strategy = strategy, test = test, study = study,
          asymptotic_p = fit$p[1L], permutation_p = fit$p[2L],
          samples = fit$samples, clusters = fit$clusters, seconds = fit$seconds,
          error = fit$error, warning = fit$warning,
          count_covariates = fit$count_covariates, zero_covariates = fit$zero_covariates,
          count_disabled = fit$count_disabled, row.names = NULL)
      }
      # Preserve the original Fisher available-study rule and report its size.
      p <- apply(values, 2L, function(z) {
        z <- z[is.finite(z)]
        if (!length(z)) return(NA_real_)
        if (any(z == 0)) return(0)
        pchisq(-2 * sum(log(z)), 2 * length(z), lower.tail = FALSE)
      })
      n.studies <- colSums(is.finite(values))
    } else {
      Y <- if (strategy == "CD") {
        if (qcat2) rare.corrected else corrected
      } else if (qcat2) rare.pooled else pooled
      Y <- Y[, intersect(colnames(Y), genera), drop = FALSE]
      fit <- run.cluster.test(Y, pooled.design, pooled.id, qcat2)
      p <- fit$p
      n.studies <- as.integer(is.finite(p)) * length(counts)
      calls[[length(calls) + 1L]] <- data.frame(lineage = family,
        strategy = strategy, test = test, study = "pooled",
        asymptotic_p = p[1L], permutation_p = p[2L],
        samples = fit$samples, clusters = fit$clusters, seconds = fit$seconds,
        error = fit$error, warning = fit$warning,
        count_covariates = fit$count_covariates, zero_covariates = fit$zero_covariates,
        count_disabled = fit$count_disabled, row.names = NULL)
    }
    records[[length(records) + 1L]] <- data.frame(lineage = family,
      strategy = strategy, test = test, calibration = c("asymptotic", "permutation"),
      method = paste0(strategy, "-", test, " (", c("asy", "perm"), ")"),
      p = unname(p), n_studies_used = unname(n.studies),
      n_perm_requested = c(0L, n.perm), n_perm_completed = NA_integer_,
      row.names = NULL)
  }
  result <- list(fingerprint = fingerprint, lineage = family,
    results = do.call(rbind, records), calls = do.call(rbind, calls),
    seconds = proc.time()[["elapsed"]] - started)
  saveRDS(result, cache.file, compress = "gzip")
  message("Completed ", family, " in ", round(result$seconds, 1L), " seconds")
  result
}

message("IBD cluster QCAT: ", length(families), " families; ", n.perm,
        " requested permutations; ", n.threads, " workers")
cluster.started <- proc.time()[["elapsed"]]
cluster.output <- parallel::mclapply(families, run.family,
  mc.cores = n.threads, mc.preschedule = FALSE,
  mc.set.seed = FALSE, mc.silent = TRUE)
if (!all(vapply(cluster.output, function(z)
  is.list(z) && is.data.frame(z$results), logical(1))))
  stop("A family worker terminated. Completed checkpoints remain available.")
cluster.results <- do.call(rbind, lapply(cluster.output, "[[", "results"))
cluster.calls <- do.call(rbind, lapply(cluster.output, "[[", "calls"))
cluster.methods <- unlist(lapply(c("P", "CD", "RD"), function(strategy)
  unlist(lapply(c("QCAT1", "QCAT2"), function(test)
    paste0(strategy, "-", test, " (", c("asy", "perm"), ")")))))
cluster.p.values <- matrix(NA_real_, length(families), length(cluster.methods),
  dimnames = list(families, cluster.methods))
for (j in seq_len(nrow(cluster.results)))
  cluster.p.values[cluster.results$lineage[j], cluster.results$method[j]] <-
    cluster.results$p[j]
cluster.q.values <- matrix(vapply(seq_len(ncol(cluster.p.values)), function(j)
  p.adjust(cluster.p.values[, j], method = "BH"), numeric(length(families))),
  nrow = length(families), dimnames = dimnames(cluster.p.values))
cluster.results$q <- cluster.q.values[cbind(
  match(cluster.results$lineage, families),
  match(cluster.results$method, cluster.methods))]
write.csv(cluster.p.values, file.path(cluster.dir, "IBD_QCAT_cluster_p_val.csv"))
write.csv(cluster.q.values, file.path(cluster.dir, "IBD_QCAT_cluster_q_val.csv"))
write.csv(cluster.results, file.path(cluster.dir, "IBD_QCAT_cluster_results.csv"),
          row.names = FALSE)
write.csv(cluster.calls,
  file.path(cluster.dir, "IBD_QCAT_cluster_call_diagnostics.csv"), row.names = FALSE)
cluster.failures <- cluster.calls[!is.finite(cluster.calls$asymptotic_p) |
  !is.finite(cluster.calls$permutation_p), , drop = FALSE]
write.csv(cluster.failures, file.path(cluster.dir, "IBD_QCAT_cluster_failures.csv"),
          row.names = FALSE)
save.result(list(configuration = configuration, p.values = cluster.p.values,
  q.values = cluster.q.values, results = cluster.results, calls = cluster.calls,
  family.results = cluster.output, preprocessing_seconds = preprocessing.seconds,
  family_wall_seconds = proc.time()[["elapsed"]] - cluster.started,
  session.info = sessionInfo()), file.path(cluster.dir, "IBD_QCAT_cluster_results.rds"))

## 8. Fresh p/q tables and five discovery sets ---------------------------------
# Asymptotic FE components, bootstrap FE/RE components, and all native comparators.
meta.methods <- c(rep(c("FE-MetaQCAT", "FE-MetaQCAT*", "FE-MetaQCAT-O"), 2L),
                  "RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O")
meta.calibrations <- c(rep("Asy", 3L), rep("Boot", 6L))
fields <- rep(c("pvalue.comp", "pvalue.pa", "pvalue.comb"), 3L)
meta.p.values <- matrix(NA_real_, length(lineages), 9L,
  dimnames = list(lineages, paste0(meta.methods, " (",
    c(rep("asy", 3L), rep("boot", 6L)), ")")))
for (j in seq_len(9L)) {
  fit <- if (j <= 6L) FE_IBD else RE_IBD
  meta.p.values[, j] <- vapply(fit$lineage.coef[lineages], function(z) {
    p <- z[[fields[j]]]
    if (is.null(p) || !meta.calibrations[j] %in% names(p)) NA_real_ else
      unname(p[meta.calibrations[j]])
  }, numeric(1))
}
meta.q.values <- matrix(vapply(seq_len(9L), function(j)
  p.adjust(meta.p.values[, j], method = "BH"), numeric(length(lineages))),
  nrow = length(lineages), dimnames = dimnames(meta.p.values))
stopifnot(isTRUE(all.equal(unname(meta.q.values[, 6L]),
  unname(FE_IBD$lineage.qval["Bootstrap-FE-MetaQCAT-O", lineages]))),
  isTRUE(all.equal(unname(meta.q.values[, 9L]),
  unname(RE_IBD$lineage.qval["Bootstrap-RE-MetaQCAT-O", lineages]))))
p.values <- cbind(meta.p.values, cluster.p.values[lineages, , drop = FALSE])
q.values <- cbind(meta.q.values, cluster.q.values[lineages, , drop = FALSE])
resampling.cols <- c(4:9, grep("\\(perm\\)$", colnames(p.values)))
write.csv(p.values, file.path(result.dir, "IBD_family_p_val.csv"))
write.csv(q.values, file.path(table.dir, "IBD_family_q_val.csv"))
write.csv(q.values[, resampling.cols, drop = FALSE],
          file.path(table.dir, "IBD_family_q_val_resampling.csv"))

lineage.rows <- list()
for (j in c(6L, 9L)) {
  component.cols <- if (j == 6L) 4:5 else 7:8
  lineage.rows[[meta.methods[j]]] <- data.frame(
    setting = "primary", method = meta.methods[j], lineage = lineages,
    p = meta.p.values[, j], q = meta.q.values[, j],
    positive_p = meta.p.values[, component.cols[1L]],
    positive_q = meta.q.values[, component.cols[1L]],
    zero_p = meta.p.values[, component.cols[2L]],
    zero_q = meta.q.values[, component.cols[2L]],
    zero_only = !is.na(meta.q.values[, j]) & meta.q.values[, j] < fdr.alpha &
      !is.na(meta.q.values[, component.cols[1L]]) &
      meta.q.values[, component.cols[1L]] >= fdr.alpha &
      !is.na(meta.q.values[, component.cols[2L]]) &
      meta.q.values[, component.cols[2L]] < fdr.alpha,
    row.names = NULL, stringsAsFactors = FALSE)
}
latest <- do.call(rbind, lineage.rows)
rownames(latest) <- NULL
write.csv(latest, file.path(bootstrap.dir, "IBD_lineage_results.csv"), row.names = FALSE)
write.csv(latest, file.path(result.dir, "IBD_lineage_results.csv"), row.names = FALSE)
fe <- lineage.rows[["FE-MetaQCAT-O"]]
re <- lineage.rows[["RE-MetaQCAT-O"]]
# Named columns replace the former fixed-position, QCAT-C-only combined table.
combined <- data.frame(lineage = lineages,
  FE_p = fe$p, FE_q = fe$q, RE_p = re$p, RE_q = re$q, check.names = FALSE)
sig <- data.frame(Lineage = lineages, check.names = FALSE)
for (strategy in c("RD", "P", "CD")) {
  columns <- paste0(strategy, "-QCAT", 1:2, " (perm)")
  sig[[strategy]] <- as.integer(rowSums(
    !is.na(cluster.q.values[, columns, drop = FALSE]) &
    cluster.q.values[, columns, drop = FALSE] < fdr.alpha) > 0L)
  for (test in 1:2) {
    column <- columns[test]
    label <- paste0(strategy, if (test == 1L) "_QCAT" else "_QCAT2")
    combined[[paste0(label, "_p")]] <- cluster.p.values[, column]
    combined[[paste0(label, "_q")]] <- cluster.q.values[, column]
  }
}
sig$FE <- as.integer(!is.na(fe$q) & fe$q < fdr.alpha)
sig$RE <- as.integer(!is.na(re$q) & re$q < fdr.alpha)
write.csv(combined, file.path(result.dir, "IBD_lineage_combine.csv"), row.names = FALSE)
write.csv(sig, file.path(result.dir, "IBD_sig_family.csv"), row.names = FALSE)
write.csv(sig, file.path(cluster.dir, "IBD_sig_family.csv"), row.names = FALSE)
discovery.summary <- data.frame(
  strategy = c("RD", "P", "CD", "FE", "RE",
               "FE_and_RE", "FE_only", "RE_only", "Any"),
  n_discoveries = c(colSums(sig[, -1L]), sum(sig$FE == 1L & sig$RE == 1L),
    sum(sig$FE == 1L & sig$RE == 0L), sum(sig$FE == 0L & sig$RE == 1L),
    sum(rowSums(sig[, -1L]) > 0L)))
write.csv(discovery.summary, file.path(bootstrap.dir, "discovery_summary.csv"),
          row.names = FALSE)

## 9. Method diagnostics, timing, configuration and session ---------------------
all.method.rows <- lapply(seq_len(ncol(p.values)), function(j) data.frame(
  lineage = lineages, method = colnames(p.values)[j],
  calibration = if (grepl("\\(asy\\)$", colnames(p.values)[j])) "Asy" else
    if (grepl("\\(boot\\)$", colnames(p.values)[j])) "Boot" else "Permutation",
  p = p.values[, j], q = q.values[, j], row.names = NULL))
write.csv(do.call(rbind, all.method.rows),
          file.path(result.dir, "IBD_all_method_results.csv"), row.names = FALSE)
run.summary <- data.frame(setting = "primary",
  method = c("FE-MetaQCAT-O", "RE-MetaQCAT-O", "Cluster-QCAT P/CD/RD"),
  samples = sum(vapply(counts, nrow, integer(1))), strata = length(counts),
  source_cohorts = length(unique(study.summary$source_cohort)),
  bootstrap_budget = c(n.boot, n.boot, 0L),
  permutation_budget = c(0L, 0L, n.perm),
  seed = seed, seconds = c(fe.seconds, re.seconds,
    proc.time()[["elapsed"]] - cluster.started),
  models_refitted = c(run.metaqcat, run.metaqcat, NA),
  failed_lineages = c(length(FE_IBD$failed.lineage), length(RE_IBD$failed.lineage),
    sum(rowSums(is.finite(cluster.p.values)) == 0L)),
  n.threads = n.threads, stringsAsFactors = FALSE)
write.csv(run.summary, file.path(bootstrap.dir, "IBD_run_summary.csv"), row.names = FALSE)
save.result(list(settings = settings, results = meta.results,
  comparators = list(configuration = configuration,
    p.values = cluster.p.values, q.values = cluster.q.values,
    results = cluster.results, calls = cluster.calls),
  p_values = p.values, q_values = q.values, significant = sig,
  study_summary = study.summary, exclusions = exclusions, run_summary = run.summary,
  elapsed_seconds = as.numeric(difftime(Sys.time(), run.started, units = "secs")),
  session_info = sessionInfo()), file.path(result.dir, "IBD_all_methods.rds"))
message("Saved IBD MetaQCAT and cluster QCAT P/CD/RD analyses under ", result.dir)
