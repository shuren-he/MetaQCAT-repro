# Make the frozen paper inputs available to the current analysis entrypoints.
# Run from the repository root. No downloads, fitting or simulations are run.
# Rscript scripts/inputs/prepare_inputs.R
# Add --check-only to validate prepared inputs without copying anything.
# Add --overwrite only to replace different files in the canonical data folders.
source("R/project_paths.R", local = TRUE)
flags <- commandArgs(trailingOnly = TRUE)
if (any(!flags %in% c("--check-only", "--overwrite")))
  stop("Supported flags: --check-only, --overwrite")
check.only <- "--check-only" %in% flags
overwrite <- "--overwrite" %in% flags

## 1. Copy the supplied, frozen inputs to the paths used by the paper scripts.
inputs <- data.frame(
  source = c(file.path("real", "CRC", "data",
                       c("count.genus.RData", "meta.RData", "tax.RData")),
             file.path("real", "IBD", "data",
                       c("count.genus.RData", "meta.RData", "tax.RData")),
             file.path("simulation", "data",
                       c("count.genus.rda", "meta.rda",
                         "mod_gdm_lst_Enterobacteriaceae.RData",
                         "mod_gdm_lst_Peptostreptococcaceae.RData"))),
  target = c(file.path("data", "CRC",
                       c("count.genus.RData", "meta.RData", "tax.RData")),
             file.path("data", "IBD",
                       c("count.genus.RData", "meta.RData", "tax.RData")),
             file.path("data", "simu",
                       c("count.genus.rda", "meta.rda",
                         "mod_gdm_lst_Enterobacteriaceae.RData",
                         "mod_gdm_lst_Peptostreptococcaceae.RData"))),
  stringsAsFactors = FALSE)
for (i in seq_len(nrow(inputs))) {
  from <- file.path(project.dir, inputs$source[i])
  to <- file.path(project.dir, inputs$target[i])
  if (!file.exists(from))
    stop("Missing bundled input: ", inputs$source[i],
         ". Obtain the complete repository, including its input files.")
  if (!check.only) {
    if (file.exists(to)) {
      same <- identical(unname(tools::md5sum(from)), unname(tools::md5sum(to)))
      if (!same && !overwrite)
        stop("Prepared input differs: ", inputs$target[i],
             ". It was not overwritten; inspect it before using --overwrite.")
    } else same <- FALSE
    if (!same) {
      dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(from, to, overwrite = overwrite))
        stop("Could not prepare: ", inputs$target[i])
    }
  }
  if (!file.exists(to)) stop("Missing prepared input: ", inputs$target[i])
}

## 2. Check the CRC counts, metadata and frozen taxonomic annotation.
crc <- new.env(parent = emptyenv())
for (name in c("count.genus.RData", "meta.RData", "tax.RData"))
  load(file.path(project.dir, "data", "CRC", name), envir = crc)
stopifnot(is.matrix(crc$count.genus),
          identical(dim(crc$count.genus), c(574L, 133L)),
          !anyNA(crc$count.genus), all(is.finite(crc$count.genus)),
          all(crc$count.genus >= 0 & crc$count.genus == floor(crc$count.genus)),
          identical(rownames(crc$count.genus), make.names(crc$meta$Sample_ID)),
          identical(colnames(crc$count.genus), unique(crc$tax$genus)),
          !anyDuplicated(crc$meta$Sample_ID), !anyDuplicated(crc$tax$OTUID),
          all(crc$meta$Group %in% c("CRC", "CTR")))

## 3. Check the IBD source tables before any analysis-level exclusions.
ibd <- new.env(parent = emptyenv())
for (name in c("count.genus.RData", "meta.RData", "tax.RData"))
  load(file.path(project.dir, "data", "IBD", name), envir = ibd)
stopifnot(is.matrix(ibd$otu), identical(dim(ibd$otu), c(249L, 4789L)),
          !anyNA(ibd$otu), all(is.finite(ibd$otu)),
          all(ibd$otu >= 0 & ibd$otu == floor(ibd$otu)),
          identical(colnames(ibd$otu), rownames(ibd$meta.data)),
          identical(rownames(ibd$otu), rownames(ibd$Tax)),
          !anyDuplicated(colnames(ibd$otu)), !anyDuplicated(rownames(ibd$otu)),
          all(c("dataset_name", "subject_accession", "sample_type", "disease",
                "gender", "race", "age") %in% colnames(ibd$meta.data)))

## 4. Check simulation templates and frozen GD parameter objects.
sim <- new.env(parent = emptyenv())
load(file.path(project.dir, "data", "simu", "count.genus.rda"), envir = sim)
load(file.path(project.dir, "data", "simu", "meta.rda"), envir = sim)
stopifnot(identical(sim$count.genus, crc$count.genus),
          identical(sim$meta, crc$meta))
study.names <- c("AUS", "CHI", "FRA", "GER", "USA")
stopifnot(identical(as.integer(table(factor(sim$meta$Country,
                                         levels = study.names))),
                    c(109L, 127L, 114L, 120L, 104L)))
for (family in c("Enterobacteriaceae", "Peptostreptococcaceae")) {
  model <- new.env(parent = emptyenv())
  load(file.path(project.dir, "data", "simu",
                 paste0("mod_gdm_lst_", family, ".RData")), envir = model)
  n.taxa <- if (family == "Enterobacteriaceae") 4L else 3L
  stopifnot(identical(names(model$mod_gdm_lst), study.names))
  for (fit in model$mod_gdm_lst) {
    stopifnot(identical(fit$col.order, seq_len(n.taxa)),
              identical(dim(fit$est.para$alpha.est), c(2L, n.taxa - 1L)),
              all(is.finite(fit$est.para$alpha.est)),
              length(fit$est.para$phi.est) == n.taxa - 1L,
              all(is.finite(fit$est.para$phi.est) & fit$est.para$phi.est > 0))
  }
  message(family, ": frozen GD parameters validated for five studies.")
}
message("Prepared and validated 10 frozen input files; no models or analyses run.")
