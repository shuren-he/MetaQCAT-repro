# Run from the repository root. No analyses are started by this script.
source("R/project_paths.R", local = TRUE)
flags <- commandArgs(trailingOnly = TRUE)
if (any(!flags %in% "--dev")) stop("Supported optional flag: --dev")
lib <- file.path(project.dir, ".library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cran <- "https://cloud.r-project.org"
core <- c("brglm2", "parallelly", "Rcpp", "RcppArmadillo", "pkgload",
          "digest", "jsonlite")
analysis <- c("MASS", "data.table", "geepack", "abind", "CompQuadForm",
              "dirmult", "ggplot2", "foreach", "doParallel", "doRNG",
              "vegan", "gplots", "rpart", "rpart.plot", "dplyr", "tidyr",
              "forestplot", "gridExtra", "ggplotify", "UpSetR", "patchwork",
              "scales")
dev <- c("testthat", "knitr", "rmarkdown", "roxygen2")
needed <- unique(c(core, analysis,
                   if ("--dev" %in% flags) dev))
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, lib = lib, repos = cran,
                   dependencies = c("Depends", "Imports", "LinkingTo"))
}
if (!requireNamespace("MMUPHin", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = lib, repos = cran)
  }
  BiocManager::install("MMUPHin", lib = lib, ask = FALSE, update = FALSE)
}
for (path in c("vendor/miLineage", "vendor/HMP")) {
  install.packages(file.path(project.dir, path), lib = lib,
                   repos = NULL, type = "source")
}
install.packages(file.path(project.dir, "packages", "MetaQCAT"), lib = lib,
                 repos = NULL, type = "source")
message("Installed sources. Primary library: ", lib)
message("No simulations, real-data analyses, or figures were run.")

