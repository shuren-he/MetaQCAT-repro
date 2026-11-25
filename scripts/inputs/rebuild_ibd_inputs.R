# Rebuild the three IBD input objects from the published MMUPHin source.
# Run from the repository root; no analysis models or results are recomputed.
# Source: Ma et al. (2022), doi:10.1186/s13059-022-02753-4.
# Upstream repository: https://github.com/biobakery/ibd_paper (MIT).
source("R/project_paths.R", local = TRUE)

## 1. Frozen public source and separate output locations -----------------------
source.commit <- "ae9d98ad170546f02af3c67797cdba971bc95e7d"
source.url <- paste0(
  "https://raw.githubusercontent.com/biobakery/ibd_paper/",
  source.commit, "/data/physeq/genera.RData")
source.dir <- Sys.getenv("INPUT_SOURCE_CACHE",
  file.path(project.dir, "data", "public", "IBD"))
output.dir <- Sys.getenv("INPUT_REBUILD_DIR",
  file.path(project.dir, "data", "rebuilt", "IBD"))
reference.dir <- file.path(project.dir, "real", "IBD", "data")
dir.create(source.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
source.file <- file.path(source.dir, "genera.RData")

## 2. Download the published, pre-filtered genus-level object -------------------
# Reuse an existing download; it must still reproduce the frozen input objects.
if (!file.exists(source.file)) {
  temporary.file <- tempfile(pattern = "genera-", tmpdir = source.dir)
  status <- utils::download.file(source.url, temporary.file, mode = "wb")
  if (status != 0L || !file.rename(temporary.file, source.file))
    stop("Could not download the public IBD source: ", source.url)
}
source.env <- new.env(parent = emptyenv())
load(source.file, envir = source.env)
if (!exists("physeq_genera", envir = source.env, inherits = FALSE))
  stop("The public source must contain physeq_genera.")
physeq <- source.env$physeq_genera
if (!identical(as.character(class(physeq)), "phyloseq"))
  stop("Unexpected class in the public IBD source.")

## 3. Extract the original counts, sample metadata and taxonomy ----------------
# The frozen source stores ordinary matrices/data.frame as phyloseq S4 slots.
# asS3() removes only the S4 flag; no phyloseq installation or fitting is needed.
# No normalization, imputation, filtering, or row/column reordering is applied.
count.slot <- attr(physeq, "otu_table")
if (!isTRUE(attr(count.slot, "taxa_are_rows")))
  stop("The frozen source must store taxa in rows.")
otu <- unclass(asS3(count.slot))
attr(otu, "taxa_are_rows") <- NULL
Tax <- unclass(asS3(attr(physeq, "tax_table")))
sample.data <- asS3(attr(physeq, "sam_data"))
class(sample.data) <- "data.frame"
attr(sample.data, ".S3Class") <- NULL
meta.data <- as.matrix(sample.data)

stopifnot(is.matrix(otu), is.matrix(Tax), is.matrix(meta.data),
          identical(dim(otu), c(249L, 4789L)),
          identical(dim(meta.data), c(4789L, 57L)),
          identical(dim(Tax), c(249L, 8L)),
          identical(colnames(otu), rownames(meta.data)),
          identical(rownames(otu), rownames(Tax)),
          !anyNA(otu), all(is.finite(otu)), all(otu >= 0),
          all(otu == floor(otu)),
          !anyDuplicated(rownames(otu)), !anyDuplicated(colnames(otu)))

## 4. Verify against the bundled analysis inputs before saving ----------------
object.names <- c("otu", "meta.data", "Tax")
file.names <- c("count.genus.RData", "meta.RData", "tax.RData")
if (!all(file.exists(file.path(reference.dir, file.names))))
  stop("The three bundled reference inputs are missing from real/IBD/data/.")
for (i in seq_along(object.names)) {
  reference <- new.env(parent = emptyenv())
  load(file.path(reference.dir, file.names[i]), envir = reference)
  if (!identical(get(object.names[i]), reference[[object.names[i]]]))
    stop("Public-source reconstruction differs from the frozen IBD object: ",
         object.names[i])
}

# Existing output files are never overwritten with different contents.
for (i in seq_along(object.names)) {
  output.file <- file.path(output.dir, file.names[i])
  if (file.exists(output.file)) {
    previous <- new.env(parent = emptyenv())
    load(output.file, envir = previous)
    if (!identical(get(object.names[i]), previous[[object.names[i]]]))
      stop("Refusing to overwrite a different rebuilt input: ", output.file)
  } else {
    save(list = object.names[i], file = output.file, compress = "xz")
  }
}
message("IBD public-source reconstruction matches all three frozen objects exactly.")
message("Saved rebuilt inputs: ", normalizePath(output.dir))
message("Primary input files and statistical results were not changed.")
