# Run analysis scripts with the repository root as the working directory.
# This file sets paths only; it does not load data or start analyses.
project.dir <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project.dir, ".metaqcat-repro-root"))) {
  stop("Set the working directory to the MetaQCAT-repro repository root first.",
       call. = FALSE)
}
script.dir <- project.dir
local.library <- file.path(project.dir, ".library")
if (dir.exists(local.library)) .libPaths(c(local.library, .libPaths()))
