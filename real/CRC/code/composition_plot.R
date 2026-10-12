# CRC real-data composition: two panels on one PDF page, matching the reference figure.
# Run: Rscript real/CRC/code/composition_plot.R
# Input: real/CRC/data/{count.genus,tax,meta}.RData
# Output: plot/real/CRC_composition.pdf (no PNGs, CSVs, or supplementary plots).
# Relative abundance is normalized within each family for each sample.
# To load without writing the PDF: options(metaqcat.plot.source.only = TRUE)

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
      if (dir.exists(file.path(candidate, "real", "CRC", "data"))) return(candidate)
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
  }
  stop("Cannot locate the MetaQCAT project; source this script by its file path.")
})


data.dir <- file.path(project.root, "real", "CRC", "data")
output.file <- file.path(project.root, "plot", "real", "CRC_composition.pdf")

make_crc_composition_panel <- function(counts, meta, tax, family, title, colors, tag) {
  genera <- unique(tax$genus[!is.na(tax$family) & tax$family == family])
  if (!length(genera) || !all(genera %in% colnames(counts))) {
    stop("Missing genera for family: ", family)
  }
  y <- counts[, genera, drop = FALSE]
  depth <- rowSums(y)
  keep <- depth > 0
  if (!any(keep)) stop("No positive counts for family: ", family)
  y <- y[keep, , drop = FALSE] / depth[keep]
  m <- meta[keep, , drop = FALSE]
  labels <- sub("^[0-9]+ ", "", genera)
  labels[grepl("^NA", genera)] <- "Unclassified"
  if (anyDuplicated(labels) || !setequal(labels, names(colors))) {
    stop("Unexpected genus labels for family: ", family)
  }
  d <- data.frame(
    sample = factor(rep(rownames(y), times = ncol(y)), levels = rownames(y)),
    Country = factor(rep(m$Country, times = ncol(y)), levels = c("AUS", "CHI", "FRA", "GER", "USA")),
    Group = factor(rep(m$Group, times = ncol(y)), levels = c("CRC", "CTR")),
    genus = factor(rep(labels, each = nrow(y)), levels = names(colors)),
    proportion = as.vector(y)
  )
  ggplot2::ggplot(d, ggplot2::aes(x = sample, y = proportion, fill = genus)) +
    ggplot2::geom_col(width = 0.88) +
    ggh4x::facet_grid2(Country ~ Group, scales = "free_x", independent = "x") +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    ggplot2::labs(title = title, tag = tag, x = NULL, y = "Relative abundance", fill = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.14, "cm"),
      panel.spacing.y = grid::unit(0.4, "cm"),
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.title.y = ggplot2::element_text(face = "bold", size = 13),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey85"),
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5),
      plot.tag = ggplot2::element_text(size = 13, face = "bold"),
      legend.position = "bottom", legend.text = ggplot2::element_text(size = 6.5),
      legend.key.width = grid::unit(0.28, "cm"), legend.key.height = grid::unit(0.3, "cm"),
      legend.spacing.x = grid::unit(0.08, "cm"),
      plot.margin = ggplot2::margin(8, 6, 8, 6)
    )
}

plot_crc_composition <- function(data.dir, output.file) {
  required <- c("ggplot2", "ggh4x", "patchwork")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  e <- new.env(parent = emptyenv())
  for (file in c("count.genus.RData", "meta.RData", "tax.RData")) {
    load(file.path(data.dir, file), envir = e)
  }
  counts <- as.matrix(e$count.genus)
  meta <- as.data.frame(e$meta, stringsAsFactors = FALSE)
  ids <- make.names(meta$Sample_ID)
  if (anyDuplicated(ids) || anyDuplicated(rownames(counts)) ||
      !setequal(rownames(counts), ids)) stop("Count and metadata sample IDs do not match.")
  meta <- meta[match(rownames(counts), ids), , drop = FALSE]
  if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      anyNA(meta[c("Country", "Group")]) ||
      !all(meta$Country %in% c("AUS", "CHI", "FRA", "GER", "USA")) ||
      !all(meta$Group %in% c("CRC", "CTR"))) stop("Invalid counts or CRC metadata.")
  # Preserve the old script's shallow-sample exclusion, using its ID rather than row 19.
  keep <- rownames(counts) != "CCIS12370844ST.4.0"
  counts <- counts[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  # Use observed counts in all countries without manually zeroing rare genera.
  pepto <- make_crc_composition_panel(counts, meta, e$tax,
    "186804 Peptostreptococcaceae", "Peptostreptococcaceae",
    c(Unclassified = "#999999", Intestinibacter = "#da6001", Peptostreptococcus = "#7570b3"), "(a)")
  entero <- make_crc_composition_panel(counts, meta, e$tax,
    "543 Enterobacteriaceae", "Enterobacteriaceae",
    c(Unclassified = "#999999", Klebsiella = "#da6001", Citrobacter = "#7570b3",
      Hafnia = "#e72b89", Enterobacter = "#67a61e"), "(b)")
  figure <- patchwork::wrap_plots(pepto, entero, nrow = 1, guides = "keep")
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(output.file, figure, width = 12, height = 10, units = "in", bg = "white")
  message("Saved CRC composition PDF: ", output.file)
  invisible(figure)
}

if (!isTRUE(getOption("metaqcat.plot.source.only", FALSE))) {
  plot_crc_composition(data.dir, output.file)
}
