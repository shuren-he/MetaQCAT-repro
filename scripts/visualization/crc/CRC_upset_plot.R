# Run from the MetaQCAT-repro repository root.
source(file.path("R", "project_paths.R"))

# Rebuild the CRC family-level UpSet plot from the significance matrix.
# Run from any directory: Rscript CRC_upset_plot.R

args <- commandArgs(trailingOnly = FALSE)
script.arg <- grep("^--file=", args, value = TRUE)
# Repository paths are defined by R/project_paths.R.

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  stop("Install UpSetR to draw this plot.", call. = FALSE)
}

sig.file <- file.path(project.dir, "res", "CRC", "sig_family.csv")
sig <- utils::read.csv(sig.file, check.names = FALSE,
                       stringsAsFactors = FALSE)
expected <- c("Lineage", "RD", "P", "CD", "FE", "RE")
if (!identical(names(sig), expected) || anyDuplicated(sig$Lineage) ||
    anyNA(sig) || !all(vapply(sig[-1L], function(x)
      is.numeric(x) && all(x %in% 0:1), logical(1)))) {
  stop("Unexpected format or values in ", sig.file, call. = FALSE)
}
sets <- lapply(sig[-1L], function(x) sig$Lineage[x == 1L])
if (any(lengths(sets) == 0L)) {
  stop("Every group must have at least one significant lineage.", call. = FALSE)
}

# Sort the sets from fewest to most discoveries, retaining the input order
# (RD, P, CD, FE, RE) to break ties. UpSetR displays its `sets` argument
# bottom-to-top when keep.order = TRUE, so reverse this display order.
set.sizes <- lengths(sets)
display.order <- names(sets)[order(set.sizes, seq_along(sets))]
upset.sets <- rev(display.order)

plot.dir <- file.path(project.dir, "plots")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
draw <- function() {
  UpSetR::upset(UpSetR::fromList(sets),
                sets = upset.sets,
                keep.order = TRUE, nintersects = 12,
                mb.ratio = c(0.5, 0.5),
                sets.x.label = "Number of Discoveries",
                order.by = c("freq", "degree"),
                decreasing = c(TRUE, FALSE),
                point.size = 3,
                text.scale = c(3, 3, 2, 3, 3, 3))
}

pdf.file <- file.path(plot.dir, "CRC_overlap_plot.pdf")
grDevices::pdf(pdf.file, width = 11, height = 8)
tryCatch(print(draw(), newpage = FALSE), finally = grDevices::dev.off())

png.file <- file.path(plot.dir, "CRC_overlap_plot.png")
grDevices::png(png.file, width = 11, height = 8,
               units = "in", res = 300, bg = "white")
tryCatch(print(draw(), newpage = FALSE), finally = grDevices::dev.off())

message("Saved CRC UpSet plot: ", pdf.file)
message("Saved CRC UpSet preview: ", png.file)
