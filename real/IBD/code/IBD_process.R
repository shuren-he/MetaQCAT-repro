rm(list=ls())

library(dirmult)
library(HMP)
library(stringr)
library(zoo)
library(qvalue)
library(grid)
library(MetaQCAT)

# Locate the project from this script (Rscript/source) or the working directory.
# This keeps paths portable if the project folder is moved or renamed.
project.root <- local({
  source.files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
  }, character(1))
  script.args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script.files <- c(rev(source.files[nzchar(source.files)]),
                    sub("^--file=", "", script.args))
  start.dirs <- unique(c(dirname(script.files), getwd()))
  for (start.dir in start.dirs) {
    if (!dir.exists(start.dir)) next
    candidate <- normalizePath(start.dir, winslash = "/", mustWork = TRUE)
    repeat {
      if (all(dir.exists(file.path(candidate, c("real/IBD/data"))))) {
        return(candidate)
      }
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
  }
  stop("Cannot locate the MetaQCAT project. Run from the project directory or source this script by its file path.")
})
data.dir <- file.path(project.root, "real", "IBD", "data")

count.env <- new.env(parent = emptyenv())
meta.env <- new.env(parent = emptyenv())
tax.env <- new.env(parent = emptyenv())
load(file.path(data.dir, "count.genus.RData"), envir = count.env)
load(file.path(data.dir, "meta.RData"), envir = meta.env)
load(file.path(data.dir, "tax.RData"), envir = tax.env)
otu <- count.env$otu
meta.data <- as.data.frame(meta.env$meta.data, stringsAsFactors = FALSE)
taxonomy <- as.matrix(tax.env$Tax[, seq_len(6L), drop = FALSE])

stopifnot(
  identical(colnames(otu), rownames(meta.data)),
  identical(rownames(otu), rownames(taxonomy))
)

# IBD_process.R intended to remove the seven shallowest samples. Select their
# sample IDs rather than using their read-depth values as row positions.
sample.depth <- colSums(otu)
low.depth.ids <- names(sort(sample.depth))[seq_len(7L)]
excluded.studies <- c("PROTECT", "LSS-PRISM", "BIDMC-FMT")
excluded.age.ids <- rownames(meta.data)[
  meta.data$dataset_name == "HMP2" &
    is.na(suppressWarnings(as.numeric(meta.data$age)))
]
excluded.samples <- data.frame(
  sample_id = c(low.depth.ids, excluded.age.ids),
  reason = c(rep("one_of_seven_lowest_depths", length(low.depth.ids)),
             rep("HMP2_age_missing", length(excluded.age.ids))),
  stringsAsFactors = FALSE
)

keep <- !rownames(meta.data) %in% excluded.samples$sample_id &
  !meta.data$dataset_name %in% excluded.studies
meta.data <- meta.data[keep, , drop = FALSE]
otu <- otu[, rownames(meta.data), drop = FALSE]
meta.data$age.numeric <- suppressWarnings(as.numeric(meta.data$age))
meta.data$analysis_study <- meta.data$dataset_name
hmp2 <- meta.data$dataset_name == "HMP2"
# Use the usual pediatric/adult terminology and retain the cohort of origin.
meta.data$analysis_study[hmp2 & meta.data$age.numeric < 18] <- "HMP2.pediatric"
meta.data$analysis_study[hmp2 & meta.data$age.numeric >= 18] <- "HMP2.adult"
stopifnot(!anyNA(meta.data$analysis_study))

meta.data$study_key <- paste(meta.data$analysis_study,
                             meta.data$sample_type, sep = ".")
metadata <- split(meta.data, meta.data$study_key, drop = TRUE)
counts <- lapply(metadata, function(m) {
  t(otu[, rownames(m), drop = FALSE])
})

# Match the original per-study prevalence rule after splitting HMP2.
prevalence.threshold <- 0.9
taxa.before.filter <- vapply(counts, ncol, integer(1))
counts <- lapply(counts, function(y) {
  y[, colMeans(y == 0) <= prevalence.threshold, drop = FALSE]
})
stopifnot(all(vapply(counts, ncol, integer(1)) > 0L))

make_design <- function(m) {
  n <- nrow(m)
  if (anyNA(m$disease)) stop("Missing disease status in ", m$study_key[1L])
  disease <- as.integer(m$disease != "control")
  X <- matrix(disease, ncol = 1L,
              dimnames = list(rownames(m), "disease"))

  if (mean(is.na(m$gender)) < 0.5) {
    # Preserve the original forward-fill rule without shortening the vector
    # when its first value is missing.
    gender <- zoo::na.locf(m$gender, na.rm = FALSE)
    gender.binary <- rep(1L, n)
    gender.binary[which(gender == "f")] <- 0L
    if (max(table(gender.binary)) / n < 0.95) {
      X <- cbind(X, gender = gender.binary)
    }
  }
  if (mean(is.na(m$race)) < 0.5) {
    race.binary <- rep(1L, n)
    race.binary[which(!is.na(m$race) & m$race != "white")] <- 0L
    if (max(table(race.binary)) / n < 0.95) {
      X <- cbind(X, race = race.binary)
    }
  }
  if (mean(is.na(m$age.numeric)) < 0.5) {
    age <- zoo::na.locf(m$age.numeric, na.rm = FALSE)
    age.binary <- rep(1L, n)
    age.binary[which(age < 18)] <- 0L
    if (max(table(age.binary)) / n < 0.95) {
      X <- cbind(X, age_adult = age.binary)
    }
  }
  X
}

design <- lapply(metadata, make_design)
subject_ids <- lapply(metadata, function(m) factor(m$subject_accession))
disease <- lapply(design, function(X) X[, 1L, drop = FALSE])
adjust <- lapply(design, function(X) X[, -1L, drop = FALSE])

for (study in names(metadata)) {
  stopifnot(
    identical(rownames(counts[[study]]), rownames(metadata[[study]])),
    identical(rownames(design[[study]]), rownames(metadata[[study]])),
    length(subject_ids[[study]]) == nrow(metadata[[study]]),
    !anyNA(subject_ids[[study]])
  )
}

summarize_age_group <- function(study, group) {
  m <- meta.data[meta.data$dataset_name == study, , drop = FALSE]
  selected <- switch(
    group,
    pediatric = !is.na(m$age.numeric) & m$age.numeric < 18,
    adult = !is.na(m$age.numeric) & m$age.numeric >= 18,
    missing = is.na(m$age.numeric)
  )
  m <- m[selected, , drop = FALSE]
  ages <- m$age.numeric[!is.na(m$age.numeric)]
  data.frame(
    study = study,
    table_label = if (study == "HMP2") "IBD4" else "IBD6",
    age_group = group,
    samples = nrow(m),
    subjects = length(unique(m$subject_accession)),
    cases = sum(m$disease != "control", na.rm = TRUE),
    controls = sum(m$disease == "control", na.rm = TRUE),
    min_age = if (length(ages)) min(ages) else NA_real_,
    max_age = if (length(ages)) max(ages) else NA_real_,
    analysis_decision = if (study == "HMP2") "split" else "retain_one_study",
    stringsAsFactors = FALSE
  )
}

age.summary <- do.call(rbind, list(
  summarize_age_group("HMP2", "pediatric"),
  summarize_age_group("HMP2", "adult"),
  summarize_age_group("MucosalIBD", "pediatric"),
  summarize_age_group("MucosalIBD", "adult")
))

study.summary <- data.frame(
  study_key = names(metadata),
  samples = vapply(metadata, nrow, integer(1)),
  subjects = vapply(metadata, function(m) length(unique(m$subject_accession)),
                    integer(1)),
  taxa_before_filter = unname(taxa.before.filter[names(metadata)]),
  taxa_after_filter = vapply(counts, ncol, integer(1)),
  row.names = NULL
)

tax = taxonomy[,5:6]
colnames(tax) <- paste0('Rank', c(1:2))

FE_Time = Sys.time()
IBD_FE <- MetaQCAT(counts, disease, 1, class.id = subject_ids, Tax = tax, n.perm = 1000,
                   meta.method = "FE-MetaQCAT-O", permute.strata = "between",
                   n.threads = 10)
FE_Time = difftime(Sys.time(), FE_Time, units = "secs")

RE_Time = Sys.time()
IBD_RE <- MetaQCAT(counts, disease, 1, class.id = subject_ids, Tax = tax, n.perm = 1000,
                   meta.method = "RE-MetaQCAT-O", permute.strata = "between",
                   n.threads = 10)
RE_Time = difftime(Sys.time(), RE_Time, units = "secs")


IBD_FE_pval <- IBD_FE$lineage.pval[2,]

# q-values for the paper's discovery plots (qvalue package).
IBD_FE_qval <- qvalue::qvalue(p = IBD_FE_pval)$qvalues


IBD_RE_pval <- IBD_RE$lineage.pval[1,]

# q-values for the paper's discovery plots (qvalue package).
IBD_RE_qval <- qvalue::qvalue(p = IBD_RE_pval)$qvalues
