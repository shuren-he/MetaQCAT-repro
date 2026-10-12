# Compare discoveries across methods without the complex real-data summary plot.
# Source this file, then supply a named list of lineage IDs (one set per method).
# Example after CRC_process.R:
#   source("real/code/upset_plot.R")
#   res <- discoveries_from_pvalues(list(FE = CRC_FE_pval, RE = CRC_RE_pval))
#   plot_discovery_upset(res, file.path(project.root, "plot", "real", "CRC_overlap.pdf"))
# IBD_FE_pval / IBD_RE_pval work the same way.
# Existing binary discovery CSVs can also be used:
#   res <- read_discovery_sets("real/CRC/outcome/sig_family.csv")
# CLI: Rscript real/code/upset_plot.R real/CRC/outcome/sig_family.csv CRC
# CLI input paths are relative to the project root; output is plot/real/CRC_overlap.pdf.
# Only PDF files are written; PNG output is not supported.

project.root <- local({
  source.files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
  }, character(1))
  args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  starts <- unique(c(dirname(c(rev(source.files[nzchar(source.files)]),
                               sub("^--file=", "", args))), getwd()))
  for (start in starts) {
    if (!dir.exists(start)) next
    candidate <- normalizePath(start, winslash = "/", mustWork = TRUE)
    repeat {
      if (dir.exists(file.path(candidate, "simulation", "data"))) return(candidate)
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
  }
  stop("Cannot locate the MetaQCAT project; source this script by its file path.")
})

# Each p-value vector must be named by lineage. BH is applied across the vector,
# rather than separately to each p-value. Use adjust.method = "none" for q-values.
discoveries_from_pvalues <- function(pvalues, alpha = 0.05, adjust.method = "BH") {
  stopifnot(is.list(pvalues), !is.null(names(pvalues)), !anyDuplicated(names(pvalues)),
            all(nzchar(names(pvalues))), length(alpha) == 1L, is.finite(alpha),
            alpha > 0, alpha < 1)
  lapply(pvalues, function(p) {
    if (!is.numeric(p) || is.null(names(p)) || anyNA(names(p)) ||
        any(!nzchar(names(p))) || anyDuplicated(names(p)) ||
        any(!is.na(p) & (!is.finite(p) | p < 0 | p > 1))) {
      stop("Provide numeric p-value vectors with unique lineage names and values in [0, 1].")
    }
    adjusted <- stats::p.adjust(p, method = adjust.method)
    names(p)[!is.na(adjusted) & adjusted < alpha]
  })
}

# First CSV column: unique lineage names; named remaining columns: 0/1 per method.
# Ignore unnamed annotation columns and rows without lineage IDs (legacy footers).
# Preserve IDs, not row positions, so overlaps reflect actual shared lineages.
read_discovery_sets <- function(file) {
  raw <- utils::read.csv(file, check.names = FALSE)
  if (ncol(raw) < 3L) stop("Expected lineage IDs followed by at least two method columns.")
  ids <- trimws(as.character(raw[[1L]]))
  keep <- !is.na(ids) & nzchar(ids)
  method.columns <- which(seq_along(raw) > 1L & nzchar(trimws(names(raw))))
  x <- raw[keep, method.columns, drop = FALSE]
  names(x) <- trimws(names(raw)[method.columns])
  if (!nrow(x) || ncol(x) < 2L || anyDuplicated(ids[keep]) ||
      anyDuplicated(names(x)) || anyNA(x) ||
      !all(vapply(x, function(column) all(column %in% c(0, 1)), logical(1)))) {
    stop("Expected a nonempty 0/1 discovery CSV with unique method names and lineage row names.")
  }
  rownames(x) <- ids[keep]
  lapply(x, function(column) rownames(x)[column == 1])
}

plot_discovery_upset <- function(res,
                                 output.file = file.path(project.root, "plot", "real", "real_data_overlap.pdf"),
                                 width = 13, height = 10) {
  if (!is.character(output.file) || length(output.file) != 1L || is.na(output.file) ||
      tolower(tools::file_ext(output.file)) != "pdf") {
    stop("output.file must be a PDF path ending in .pdf.")
  }
  if (!requireNamespace("UpSetR", quietly = TRUE)) stop("Install UpSetR to draw this plot.")
  if (!is.list(res) || is.data.frame(res) || length(res) < 2L ||
      is.null(names(res)) || anyNA(names(res)) || any(!nzchar(names(res))) ||
      anyDuplicated(names(res))) {
    stop("res must be a named list of lineage vectors for at least two methods.")
  }
  if (!all(vapply(res, function(x) is.atomic(x) || is.null(x), logical(1)))) {
    stop("Each discovery set must be a vector of lineage IDs.")
  }
  res <- lapply(res, function(x) unique(as.character(x[!is.na(x)])))
  if (any(vapply(res, function(x) any(!nzchar(x)), logical(1)))) stop("Empty lineage IDs are not allowed.")
  if (!length(unique(unlist(res, use.names = FALSE)))) {
    message("No discoveries in any method; no UpSet plot was created.")
    return(invisible(NULL))
  }
  # UpSetR cannot display zero-size sets; retain their names in a message.
  empty <- lengths(res) == 0L
  if (any(empty)) message("Methods with no discoveries: ", paste(names(res)[empty], collapse = ", "))
  nonempty <- res[!empty]
  if (length(nonempty) < 2L) {
    message("At least two methods need discoveries to draw an UpSet comparison.")
    return(invisible(NULL))
  }
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(output.file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot <- UpSetR::upset(
    UpSetR::fromList(nonempty), sets = names(nonempty), nsets = length(nonempty),
    nintersects = 12, mb.ratio = c(0.5, 0.5),
    sets.x.label = "Number of Discoveries", order.by = c("freq", "degree"),
    decreasing = c(TRUE, FALSE), point.size = 3,
    text.scale = c(3, 3, 2, 3, 3, 3)
  )
  print(plot)
  message("Saved UpSet PDF: ", output.file)
  invisible(plot)
}

# Run the CLI only when this file itself was invoked by Rscript, not when sourced.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    message("Supply a project-relative discovery CSV and optional label, or source this file and call plot_discovery_upset(res).")
  } else {
    label <- if (length(args) >= 2L) args[2L] else "real_data"
    label <- gsub("[^A-Za-z0-9_-]", "_", label)
    input <- file.path(project.root, args[1L])
    plot_discovery_upset(read_discovery_sets(input),
                         file.path(project.root, "plot", "real", paste0(label, "_overlap.pdf")))
  }
}
