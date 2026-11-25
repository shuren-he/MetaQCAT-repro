# Rebuild the CRC count/metadata inputs from the versioned public profiles.
# Run from the repository root; ordinary sections can also be run in RStudio.
# The supplied taxonomy is a frozen annotation mapping, not inferred from counts.
# No association tests, rarefaction, or generalized-Dirichlet fits are run.
source("R/project_paths.R", local = TRUE)
if (!requireNamespace("data.table", quietly = TRUE))
  stop("Install data.table before rebuilding the public CRC inputs.")

## 1. Locations: never overwrite the supplied or canonical analysis inputs ----
frozen.dir <- file.path(project.dir, "real", "CRC", "data")
source.dir <- Sys.getenv("INPUT_SOURCE_CACHE",
  file.path(project.dir, "data", "public", "CRC"))
output.dir <- Sys.getenv("INPUT_REBUILD_DIR",
  file.path(project.dir, "data", "rebuilt", "CRC"))
output.names <- c("count.genus.RData", "meta.RData", "tax.RData")
if (any(file.exists(file.path(output.dir, output.names))))
  stop("Rebuilt input files already exist; choose a fresh INPUT_REBUILD_DIR.")
if (!all(file.exists(file.path(frozen.dir, output.names))))
  stop("Restore the three supplied frozen inputs under real/CRC/data/.")
frozen <- new.env(parent = emptyenv())
for (name in output.names) load(file.path(frozen.dir, name), envir = frozen)

## 2. Download the two frozen public source tables, or reuse a verified cache --
# Wirbel et al., Zenodo version 1.0, DOI 10.5281/zenodo.3517209.
source.names <- c("species_profiles_g2_l75_motus2.0.0.tsv", "meta_all.tsv")
source.urls <- paste0("https://zenodo.org/api/records/3517209/files/",
                      source.names, "/content")
source.md5 <- c("7c754ad66c39d4f13da4d0df1549bd76",
                "da18b10fdabae6308329e80b73991f84")
dir.create(source.dir, recursive = TRUE, showWarnings = FALSE)
source.files <- file.path(source.dir, source.names)
for (i in seq_along(source.files)) {
  if (!file.exists(source.files[i]))
    utils::download.file(source.urls[i], source.files[i], mode = "wb")
}
stopifnot(identical(unname(tools::md5sum(source.files)), source.md5))

## 3. Select the original five studies before reproducing the feature screen ---
upstream.meta <- data.table::fread(source.files[2L], data.table = FALSE,
  check.names = FALSE, showProgress = FALSE)
source.studies <- c("AT-CRC", "CN-CRC", "FR-CRC", "DE-CRC", "US-CRC")
source.meta <- upstream.meta[upstream.meta$Study %in% source.studies &
  upstream.meta$Group %in% c("CRC", "CTR"), , drop = FALSE]
stopifnot(nrow(source.meta) == 575L, !anyDuplicated(source.meta$Sample_ID))
# The TSV header has sample names but no heading for its feature-name column.
header <- strsplit(readLines(source.files[1L], n = 1L, warn = FALSE),
                   "\t", fixed = TRUE)[[1L]]
stopifnot(all(source.meta$Sample_ID %in% header))
profiles <- data.table::fread(source.files[1L], sep = "\t", skip = 1L,
  header = FALSE, select = c(1L, 1L + match(source.meta$Sample_ID, header)),
  col.names = c("feature", source.meta$Sample_ID), data.table = FALSE,
  check.names = FALSE, showProgress = FALSE, nThread = 2L)
otu.id <- sub("^.*\\[([^][]+)\\]$", "\\1", profiles$feature)
profiles <- as.matrix(profiles[, source.meta$Sample_ID, drop = FALSE])
rownames(profiles) <- otu.id
stopifnot(!anyDuplicated(otu.id), all(is.finite(profiles)), all(profiles >= 0))
# Relative abundance is used only for the upstream species screen. The
# denominator includes the unassigned (-1) row, following prepare_data.R.
relative <- sweep(profiles, 2L, colSums(profiles), "/")
max.by.study <- vapply(source.studies, function(study) {
  apply(relative[, source.meta$Sample_ID[source.meta$Study == study],
                 drop = FALSE], 1L, max)
}, numeric(nrow(profiles)))
selected.ids <- otu.id[rowSums(max.by.study >= 0.001) >= 3L & otu.id != "-1"]
tax <- frozen$tax
stopifnot(length(selected.ids) == 849L,
          identical(selected.ids, as.character(tax$OTUID)))

## 4. Restore the frozen sample universe, then sum original counts by genus ----
# The supplied analysis input excludes ERR1018294. Its historical exclusion
# reason was not recorded; it is NOT the later depth <1,000 exclusion.
# Keep this explicit frozen choice so existing analyses are not changed.
excluded.source.sample <- "ERR1018294"
stopifnot(excluded.source.sample %in% source.meta$Sample_ID)
meta <- source.meta[source.meta$Sample_ID != excluded.source.sample, , drop = FALSE]
rownames(meta) <- NULL
# Match numeric storage of the supplied metadata; no values are imputed.
numeric.fields <- c("Age", "BMI", "Diabetes", "Vegetarian", "Smoking", "Library_Size")
meta[numeric.fields] <- lapply(meta[numeric.fields], as.numeric)
meta$block <- ifelse(meta$Study == "CN-CRC",
  paste(meta$Study, meta$Sampling_rel_to_colonoscopy, sep = "_"), meta$Study)
stopifnot(nrow(meta) == 574L, identical(meta$Sample_ID, frozen$meta$Sample_ID),
          identical(names(meta), names(frozen$meta)),
          all(vapply(names(meta), function(name)
            identical(meta[[name]], frozen$meta[[name]]), logical(1))))
selected <- profiles[match(selected.ids, otu.id), meta$Sample_ID, drop = FALSE]
# These are count-scale marker profiles. Do not multiply relative profiles by
# a library size, round again, or rarefy while reconstructing the input table.
count.genus <- t(rowsum(selected, group = tax$genus, reorder = FALSE))
rownames(count.genus) <- make.names(meta$Sample_ID)
stopifnot(identical(count.genus, frozen$count.genus),
          identical(dim(count.genus), c(574L, 133L)))
# Metadata values/column types match; the rebuilt object is a plain data.frame
# rather than requiring the historical readr/tibble printing attributes.

## 5. Save separate rebuilt inputs; analysis data and results stay untouched ----
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
save(count.genus, file = file.path(output.dir, "count.genus.RData"), compress = "xz")
save(meta, file = file.path(output.dir, "meta.RData"), compress = "xz")
save(tax, file = file.path(output.dir, "tax.RData"), compress = "xz")
writeLines(c(
  "CRC source-to-input reconstruction",
  "Public source: https://doi.org/10.5281/zenodo.3517209 (version 1.0)",
  source.urls,
  "Feature screen: max relative abundance >=0.001 in >=3 of 5 studies; exclude -1.",
  "Screen uses 575 source samples; frozen analysis universe excludes ERR1018294.",
  "Taxonomy: supplied frozen real/CRC/data/tax.RData annotation mapping.",
  "849 mOTUs aggregated to 133 genera; 574 samples; count matrix identical.",
  "All metadata values and column types identical; rebuilt metadata is a plain data.frame.",
  "No model fitting, sample-depth filtering, or association tests were performed."
), file.path(output.dir, "reconstruction.txt"))
message("Rebuilt and verified CRC inputs: ", normalizePath(output.dir))
