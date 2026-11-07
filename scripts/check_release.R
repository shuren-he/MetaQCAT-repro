# Static checks only: no models, simulations, downloads, or installations.
source("R/project_paths.R", local = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--write-checksums")) stop("Unknown check flag.")

files <- sort(list.files(project.dir, recursive = TRUE, all.files = TRUE,
                         full.names = FALSE))
files <- files[!file.info(file.path(project.dir, files))$isdir]
files <- files[!grepl("^([.]git|[.]library|tmp|output|outputs)/", files)]
r.files <- files[grepl("[.]R$", files)]
errors <- character()
for (path in r.files) {
  error <- tryCatch({ parse(file.path(project.dir, path)); NULL },
                    error = function(e) conditionMessage(e))
  if (!is.null(error)) errors <- c(errors, paste(path, error, sep = ": "))
}
if (length(errors)) stop(paste(errors, collapse = "\n"))

artifacts <- files[grepl("[.](o|so|dll|dylib|a)$|(^|/)[.]DS_Store$|(^|/)[.]Rhistory$",
                         files)]
if (length(artifacts)) stop("Build/session artifacts found: ",
                            paste(artifacts, collapse = ", "))
external <- files[grepl("^(data|res|outcome|simu_res|type1_error)/", files)]
if (length(external)) stop("External analysis inputs/outputs found: ",
                           paste(external, collapse = ", "))
old.workflows <- files[grepl("(^|/)legacy(/|_)|^packages/[^/]+/", files) &
  !grepl("^packages/MetaQCAT/", files)]
if (length(old.workflows)) stop("Non-paper or historical workflow found: ",
                                paste(old.workflows, collapse = ", "))
entrypoints <- files[grepl("^scripts/real_data/.+[.]R$", files)]
expected <- c("scripts/real_data/crc/CRC_all_methods.R",
              "scripts/real_data/ibd/IBD_all_methods.R")
if (!setequal(entrypoints, expected))
  stop("Keep exactly one CRC and one IBD real-data entrypoint.")
binary.data <- files[grepl("[.](rda|RData|rds|Rdata)$", files)]
outside.packages <- binary.data[!grepl("^(packages|vendor)/", binary.data)]
if (length(outside.packages)) stop("Non-package data found: ",
                                   paste(outside.packages, collapse = ", "))

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Install digest to verify SHA256SUMS (or run install.R).")
}
hashed <- setdiff(files, "SHA256SUMS")
hash <- vapply(file.path(project.dir, hashed), digest::digest, character(1),
                algo = "sha256", file = TRUE, serialize = FALSE)
checksum.lines <- paste(hash, hashed, sep = "  ")
manifest <- file.path(project.dir, "SHA256SUMS")
if ("--write-checksums" %in% args) {
  writeLines(checksum.lines, manifest, useBytes = TRUE)
} else {
  if (!file.exists(manifest)) stop("Missing SHA256SUMS.")
  recorded <- readLines(manifest, warn = FALSE)
  if (!identical(recorded, unname(checksum.lines))) {
    stop("The release differs from SHA256SUMS. Review edits before refreshing it.")
  }
}
cat("R files parsed:", length(r.files), "\n")
cat("Files checksummed:", length(hashed), "\n")
cat("Package example-data files:", length(binary.data), "\n")
cat("No external analysis data/results or compiled/session artifacts found.\n")

