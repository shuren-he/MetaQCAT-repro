# Stage user-generated summaries for the current power visualization.
# Run from the repository root after completing the simulation workflows.
# Existing res/ files are preserved unless --overwrite is supplied explicitly.
source("R/project_paths.R", local = TRUE)

families <- c("Peptostreptococcaceae", "Enterobacteriaceae")
sources <- character()
for (family in families) {
  sources <- c(sources, file.path(project.dir, "outcome", "simulation", "power",
    family, paste0(family, "_power_summary.csv")))
}
missing <- sources[!file.exists(sources)]
if (length(missing)) {
  stop("Run the relevant simulations or restore their summaries first:\n",
       paste(missing, collapse = "\n"), call. = FALSE)
}
destinations <- file.path(project.dir, "res", basename(sources))
overwrite <- "--overwrite" %in% commandArgs(trailingOnly = TRUE)
if (!overwrite && any(file.exists(destinations))) {
  stop("Some res/ summaries already exist. Review them before rerunning with --overwrite.",
       call. = FALSE)
}
dir.create(file.path(project.dir, "res"), showWarnings = FALSE, recursive = TRUE)
copied <- file.copy(sources, destinations, overwrite = overwrite)
if (!all(copied)) stop("Could not stage every summary.", call. = FALSE)
message("Staged two all-method power summaries in res/. No statistical analyses were run.")
