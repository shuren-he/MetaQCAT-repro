# Run from the MetaQCAT-repro repository root.
source(file.path("R", "project_paths.R"))

# Rebuild the IBD UpSet plot from bootstrap FE/RE and existing QCAT-C results.
# Run from any directory: Rscript IBD_upset_plot_bootstrap.R

args <- commandArgs(trailingOnly = FALSE)
script.arg <- grep("^--file=", args, value = TRUE)
# Repository paths are defined by R/project_paths.R.

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  stop("Install UpSetR to draw this plot.", call. = FALSE)
}

result.file <- file.path(project.dir, "res", "IBD", "IBD_lineage_combine.csv")
result <- utils::read.csv(result.file, stringsAsFactors = FALSE)
if (ncol(result) != 15L ||
    !identical(names(result)[c(1L, 3L, 5L, 8L, 10L, 14L, 15L)],
               c("lineage", "method", "q", "method.1", "q.1", "p.2", "q.2"))) {
  stop("Unexpected column layout in ", result.file, call. = FALSE)
}

# Ignore only rows with
# neither a lineage name nor test results; never turn a missing q-value into 0.
named <- !is.na(result$lineage) & nzchar(trimws(result$lineage))
if (any(!named & !is.na(result[[5L]])) ||
    any(!named & !is.na(result[[10L]])) ||
    any(!named & !is.na(result[[15L]]))) {
  stop("An unlabelled row contains a q-value in ", result.file, call. = FALSE)
}
result <- result[named, , drop = FALSE]
if (nrow(result) != 45L || anyDuplicated(result$lineage) ||
    !all(result[[3L]] == "FE-MetaQCAT-O") ||
    !all(result[[8L]] == "RE-MetaQCAT-O") ||
    !all(result[[13L]] == "QCAT-C")) {
  stop("Unexpected lineages or method labels in ", result.file, call. = FALSE)
}

q.values <- result[c(15L, 5L, 10L)]
names(q.values) <- c("QCAT-C", "FE", "RE")
if (!all(vapply(q.values, is.numeric, logical(1))) ||
    any(unlist(q.values, use.names = FALSE) < 0 |
        unlist(q.values, use.names = FALSE) > 1, na.rm = TRUE)) {
  stop("Invalid q-values in ", result.file, call. = FALSE)
}

sig <- data.frame(Lineage = result$lineage, check.names = FALSE)
for (method in names(q.values)) {
  q <- q.values[[method]]
  sig[[method]] <- as.integer(!is.na(q) & q < 0.05)
}
sig.file <- file.path(project.dir, "res", "IBD", "IBD_sig_family.csv")
utils::write.csv(sig, sig.file, row.names = FALSE)

sets <- lapply(sig[-1L], function(x) sig$Lineage[x == 1L])
if (any(lengths(sets) == 0L)) {
  stop("Every group must have at least one significant lineage.", call. = FALSE)
}

# Display top-to-bottom by increasing discovery count. UpSetR renders its
# `sets` argument bottom-to-top when keep.order = TRUE.
display.order <- names(sets)[order(lengths(sets), seq_along(sets))]
upset.sets <- rev(display.order)
draw <- function() {
  UpSetR::upset(UpSetR::fromList(sets),
                sets = upset.sets, keep.order = TRUE,
                nintersects = 12, mb.ratio = c(0.5, 0.5),
                sets.x.label = "Number of Discoveries",
                order.by = c("freq", "degree"),
                decreasing = c(TRUE, FALSE), point.size = 3,
                text.scale = c(3, 3, 2, 3, 3, 3))
}

plot.dir <- file.path(project.dir, "outcome", "IBD")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
pdf.file <- file.path(plot.dir, "IBD_overlap_plot.pdf")
grDevices::pdf(pdf.file, width = 11, height = 8)
tryCatch(print(draw(), newpage = FALSE), finally = grDevices::dev.off())

png.file <- file.path(plot.dir, "IBD_overlap_plot.png")
grDevices::png(png.file, width = 11, height = 8,
               units = "in", res = 300, bg = "white")
tryCatch(print(draw(), newpage = FALSE), finally = grDevices::dev.off())

message("Saved IBD significance matrix: ", sig.file)
message("Saved IBD UpSet plot: ", pdf.file)
message("Saved IBD UpSet preview: ", png.file)
