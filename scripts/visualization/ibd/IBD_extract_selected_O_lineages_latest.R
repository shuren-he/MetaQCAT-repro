# Run from the MetaQCAT-repro repository root.
source(file.path("R", "project_paths.R"))

# Extract taxon-level summaries for the current nonoverlapping IBD analysis.
# The two families are the ones used in the paper's selected-lineage figures;
# their current significance is recorded in selection.csv, not assumed.

script.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
# Repository paths are defined by R/project_paths.R.

families <- c("f__Ruminococcaceae", "f__Peptostreptococcaceae")
requested.families <- commandArgs(trailingOnly = TRUE)
if (length(requested.families)) families <- requested.families
if (anyDuplicated(families)) stop("Selected family names must be unique.")
default.result.dir <- file.path(project.dir, "outcome", "IBD",
                                "selected_O_lineages")
selected.dir <- Sys.getenv("IBD_SELECTED_DIR", unset = "")
result.dir <- if (!nzchar(selected.dir)) default.result.dir else if (
  startsWith(selected.dir, "/")) selected.dir else
    file.path(project.dir, selected.dir)
data.file <- file.path(project.dir, "res", "IBD", "wild_score_bootstrap",
                       "IBD_processed_bootstrap.RData")
combined.file <- file.path(project.dir, "res", "IBD",
                           "IBD_lineage_combine.csv")
fe.file <- file.path(project.dir, "res", "IBD", "wild_score_bootstrap",
                     "IBD_primary_FE_MetaQCAT_O.rds")
re.file <- file.path(project.dir, "res", "IBD", "wild_score_bootstrap",
                     "IBD_primary_RE_MetaQCAT_O.rds")

dat <- new.env(parent = emptyenv())
load(data.file, envir = dat)
study.names <- names(dat$counts)
if (!identical(study.names, names(dat$design)) ||
    !identical(study.names, names(dat$subject_ids)) ||
    !all(vapply(study.names, function(study) {
      identical(rownames(dat$counts[[study]]),
                rownames(dat$design[[study]])) &&
        nrow(dat$counts[[study]]) == length(dat$subject_ids[[study]]) &&
        !anyNA(dat$design[[study]])
    }, logical(1)))) {
  stop("Current IBD counts, adjusted design and subject IDs are not aligned.")
}
if (!all(vapply(dat$design, function(x) "disease" %in% colnames(x),
                logical(1))) ||
    !all(vapply(dat$design, function(x) colnames(x)[[1L]] == "disease",
                logical(1)))) {
  stop("The tested disease variable must be the first adjusted-design column.")
}

combined <- utils::read.csv(combined.file, check.names = TRUE)
if (ncol(combined) != 15L ||
    anyDuplicated(combined$lineage[nzchar(combined$lineage)]) ||
    !all(families %in% combined$lineage) ||
    !all(combined$method[nzchar(combined$lineage)] == "FE-MetaQCAT-O") ||
    !all(combined$method.1[nzchar(combined$lineage)] == "RE-MetaQCAT-O")) {
  stop("Unexpected latest IBD combined-lineage table format.")
}
fe <- readRDS(fe.file)
re <- readRDS(re.file)
fe.row <- "Bootstrap-FE-MetaQCAT-O"
re.row <- "Bootstrap-RE-MetaQCAT-O"
if (!all(families %in% colnames(fe$lineage.pval)) ||
    !all(families %in% colnames(re$lineage.pval)) ||
    !fe.row %in% rownames(fe$lineage.pval) ||
    !re.row %in% rownames(re$lineage.pval)) {
  stop("Selected families are missing from the current IBD FE/RE results.")
}
selected <- combined[match(families, combined$lineage), ]
close.enough <- function(x, y) all(is.finite(x) & is.finite(y) &
                                    abs(x - y) < 1e-7)
if (!close.enough(selected$p, fe$lineage.pval[fe.row, families]) ||
    !close.enough(selected$q, fe$lineage.qval[fe.row, families]) ||
    !close.enough(selected$p.1, re$lineage.pval[re.row, families]) ||
    !close.enough(selected$q.1, re$lineage.qval[re.row, families])) {
  stop("Combined IBD p/q values differ from the current FE/RE RDS results.")
}

meta.sub <- new.env(parent = globalenv())
sys.source(file.path(project.dir, "code", "MetaSub.R"), envir = meta.sub)
if (!exists("MetaSub_latest", envir = meta.sub, inherits = FALSE)) {
  stop("MetaSub_latest() is missing from code/MetaSub.R.")
}
dir.create(result.dir, recursive = TRUE, showWarnings = FALSE)

effect.tables <- vector("list", length(families))
for (i in seq_along(families)) {
  family <- families[[i]]
  message("Extracting latest IBD summary: ", family)
  stem <- sub("^f__", "", family)
  reuse.file <- file.path(default.result.dir,
                          paste0(stem, "_latest_summary.rds"))
  reuse <- FALSE
  if (family == "f__Coriobacteriaceae" &&
      !identical(result.dir, default.result.dir) && file.exists(reuse.file)) {
    prior <- readRDS(reuse.file)
    reuse <- identical(prior$lineage, family) &&
      identical(prior$min.prev, 0.1) &&
      identical(prior$min.depth, 0) &&
      identical(prior$seed, 1234L) &&
      identical(prior$analysis_source, basename(data.file)) &&
      isTRUE(prior$analysis_samples ==
               sum(vapply(dat$counts, nrow, integer(1)))) &&
      identical(prior$positive$reference.taxon,
                fe$lineage.coef[[family]]$reference.taxon) &&
      identical(names(prior$counts), study.names) &&
      identical(names(prior$rarefied_counts), study.names) &&
      all(vapply(study.names, function(study) {
        identical(rownames(prior$counts[[study]]),
                  rownames(dat$counts[[study]])) &&
          identical(rownames(prior$rarefied_counts[[study]]),
                    rownames(dat$counts[[study]]))
      }, logical(1)))
  }
  if (reuse) {
    message("Reusing validated summary: ", reuse.file)
    summary <- prior
  } else {
    summary <- meta.sub$MetaSub_latest(
      counts = dat$counts, X = dat$design, X.index = 1L,
      taxonomy = dat$taxonomy, lineage = family,
      cluster.id = dat$subject_ids, min.depth = 0,
      min.prev = 0.1, seed = 1234L
    )
  }
  if (!identical(names(summary$counts), study.names) ||
      !identical(names(summary$rarefied_counts), study.names) ||
      !all(vapply(study.names, function(study) {
        identical(rownames(summary$counts[[study]]),
                  rownames(dat$counts[[study]])) &&
          identical(rownames(summary$rarefied_counts[[study]]),
                    rownames(dat$counts[[study]]))
      }, logical(1)))) {
    stop("MetaSub sample/study alignment failed for ", family, ".")
  }
  current.reference <- fe$lineage.coef[[family]]$reference.taxon
  if (!identical(summary$positive$reference.taxon, current.reference)) {
    stop("MetaSub reference taxon differs from the current FE result for ",
         family, ": ", summary$positive$reference.taxon, " vs ",
         current.reference)
  }
  summary$analysis_source <- basename(data.file)
  summary$analysis_samples <- sum(vapply(dat$counts, nrow, integer(1)))
  summary$analysis_design <- "disease + observed adjustment covariates"
  saveRDS(summary, file.path(result.dir, paste0(stem,
                                               "_latest_summary.rds")),
          compress = "xz")

  effect.tables[[i]] <- do.call(rbind, lapply(c("positive", "zero"),
    function(part) {
      out <- summary[[part]]$estimates
      data.frame(lineage = family, part = part, out,
                 row.names = NULL, check.names = FALSE)
    }))
}

taxon.effects <- do.call(rbind, effect.tables)
rownames(taxon.effects) <- NULL
if (anyNA(taxon.effects) ||
    anyDuplicated(taxon.effects[c("lineage", "part", "study", "taxon")])) {
  stop("The latest IBD taxon-effect table has missing or duplicate estimates.")
}

taxa <- dat$taxonomy
observed.taxa <- unique(unlist(lapply(dat$counts, colnames), use.names = FALSE))
genus.count <- vapply(families, function(family) {
  length(unique(taxa[rownames(taxa) %in% observed.taxa &
                       taxa[, 5L] == family, 6L]))
}, integer(1))
selection <- data.frame(
  lineage = families,
  genus_level_taxa = unname(genus.count),
  fe_resampling_p = selected$p,
  fe_resampling_q = selected$q,
  re_resampling_p = selected$p.1,
  re_resampling_q = selected$q.1,
  latest_significance = ifelse(selected$q < 0.05 & selected$q.1 < 0.05,
                               "FE and RE", ifelse(selected$q < 0.05,
                                                   "FE only", ifelse(
                                                     selected$q.1 < 0.05,
                                                     "RE only", "neither"))),
  sample_count = sum(vapply(dat$counts, nrow, integer(1))),
  source_data = basename(data.file),
  min_prev = 0.1,
  rarefaction_seed = 1234L,
  row.names = NULL,
  check.names = FALSE
)
if (!identical(result.dir, default.result.dir)) {
  figure.roles <- c(
    "f__Ruminococcaceae" = "positive and zero signals",
    "f__Coriobacteriaceae" = "broad genus coverage",
    "f__Pasteurellaceae" = "genus-level observed-zero illustration",
    "f__Actinomycetaceae" = "sparse genera; supplementary illustration"
  )
  selection$figure_role <- unname(figure.roles[selection$lineage])
}
utils::write.csv(selection, file.path(result.dir, "selection.csv"),
                 row.names = FALSE)
utils::write.csv(taxon.effects,
                 file.path(result.dir, "taxon_effects.csv"), row.names = FALSE)
message("Saved latest IBD summaries and taxon effects in: ", result.dir)
