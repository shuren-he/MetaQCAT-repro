# Static checks only: no models, simulations, downloads, or installations.
source("R/project_paths.R", local = TRUE)

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
# These ten frozen source inputs are required release assets.
# This is not a blanket exception for output RDS files or arbitrary data.
allowed.inputs <- c(
  "real/CRC/data/count.genus.RData", "real/CRC/data/meta.RData",
  "real/CRC/data/tax.RData", "real/IBD/data/count.genus.RData",
  "real/IBD/data/meta.RData", "real/IBD/data/tax.RData",
  "simulation/data/count.genus.rda", "simulation/data/meta.rda",
  "simulation/data/mod_gdm_lst_Enterobacteriaceae.RData",
  "simulation/data/mod_gdm_lst_Peptostreptococcaceae.RData")
missing.inputs <- allowed.inputs[!file.exists(file.path(project.dir, allowed.inputs))]
if (length(missing.inputs))
  stop("Missing frozen release inputs: ", paste(missing.inputs, collapse = ", "))
binary.data <- files[grepl("[.](rda|RData|rds|Rdata)$", files)]
outside.packages <- binary.data[!grepl("^(packages|vendor)/", binary.data) &
                                  !binary.data %in% allowed.inputs]
if (length(outside.packages)) stop("Non-package data found: ",
                                   paste(outside.packages, collapse = ", "))

cat("R files parsed:", length(r.files), "\n")
cat("Package/vendor data files:",
    sum(grepl("^(packages|vendor)/", binary.data)), "\n")
cat("Bundled frozen analysis inputs:",
    sum(binary.data %in% allowed.inputs), "\n")
cat("No unlisted data, analysis results or compiled/session artifacts found.\n")

