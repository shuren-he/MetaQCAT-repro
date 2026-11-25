# Correlated-sample power: all cluster comparator permutations + MetaQCAT bootstrap.
# Run from the repository root:
# Rscript scripts/visualization/simulation/simu_cluster_power_bar_plot.R
# Same layout/style as simu_power_bar_plot.R; independent figures stay unchanged.
# QCAT1 = QCAT.Cluster; QCAT2 = QCAT_GEE.Cluster. No DM or ZIGDM.
# Plot the supplied power = rejections / all attempted replicates, not valid-only power.

required <- c("ggplot2", "patchwork", "scales", "gtable")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

source("R/project_paths.R", local = TRUE)
project_dir <- project.dir

input_dir <- file.path(project_dir, "outcome", "simulation", "power_cluster")
plot_dir <- file.path(project_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Use the method order and group spacing from bar_plot.R.
method_positions <- c(
  "P-QCAT1 (perm)" = 1,
  "P-QCAT2 (perm)" = 2,
  "CD-QCAT1 (perm)" = 4,
  "CD-QCAT2 (perm)" = 5,
  "RD-QCAT1 (perm)" = 7,
  "RD-QCAT2 (perm)" = 8,
  "FE-MetaQCAT (boot)" = 10,
  "FE-MetaQCAT* (boot)" = 11,
  "FE-MetaQCAT-O (boot)" = 12,
  "RE-MetaQCAT (boot)" = 14,
  "RE-MetaQCAT* (boot)" = 15,
  "RE-MetaQCAT-O (boot)" = 16
)

taxa <- c("Peptostreptococcaceae", "Enterobacteriaceae")
read_power <- function(taxon) {
  file <- file.path(input_dir, taxon, paste0(taxon,
    "_correlated_power_resampling_summary.csv"))
  if (!file.exists(file)) stop("Missing correlated summary: ", file, call. = FALSE)
  x <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("lineage", "design", "scenario", "signal.method", "signal.size",
              "method", "calibration", "power", "n_rejections", "n_replicates",
              "n_valid", "failures", "n_boot_requested", "n_perm_requested")
  if (!all(needed %in% names(x)) || anyDuplicated(names(x)))
    stop("Unexpected summary columns in: ", file, call. = FALSE)
  if (anyNA(x[needed[c(1:7)]]) || any(x$lineage != taxon) ||
      any(x$design != "correlated"))
    stop("Unexpected lineage, design or missing setting labels in: ", file,
         call. = FALSE)

  # The two CSVs already contain resampling results only.
  is_meta <- grepl("MetaQCAT", x$method, fixed = TRUE)
  if (anyNA(x$calibration) ||
      any(x$calibration[is_meta] != "wild-score bootstrap") ||
      any(x$calibration[!is_meta] != "between-cluster permutation"))
    stop("Unexpected calibration in: ", file, call. = FALSE)
  x$source_method <- x$method
  x$method <- sub(" (perm)", "", x$method, fixed = TRUE)
  x$method <- sub("QCAT_GEE.Cluster", "QCAT2", x$method, fixed = TRUE)
  x$method <- sub("QCAT.Cluster", "QCAT1", x$method, fixed = TRUE)
  x$method <- paste0(x$method, ifelse(is_meta, " (boot)", " (perm)"))

  expected <- expand.grid(
    method = names(method_positions),
    signal.method = c("proportion", "count"),
    scenario = c("homogeneous", "one.flipped.study"),
    stringsAsFactors = FALSE
  )
  key <- function(d) paste(d$method, d$signal.method, d$scenario, sep = "\r")
  if (nrow(x) != nrow(expected) || anyDuplicated(key(x)) ||
      !setequal(key(x), key(expected)))
    stop("Expected 12 methods in each of four correlated settings in: ", file,
         call. = FALSE)
  x <- x[match(key(expected), key(x)), , drop = FALSE]

  # Preserve the supplied all-attempted denominator for EVERY method.
  counts <- x[c("n_replicates", "n_valid", "n_rejections", "failures",
                 "n_boot_requested", "n_perm_requested")]
  if (anyNA(counts) || any(!is.finite(as.matrix(counts))) ||
      any(as.matrix(counts) != floor(as.matrix(counts))) ||
      any(!is.finite(x$power) | x$power < 0 | x$power > 1 |
          x$n_replicates < 1L | x$n_valid < 0L |
          x$n_valid > x$n_replicates | x$n_rejections < 0L |
          x$n_rejections > x$n_valid |
          x$failures != x$n_replicates - x$n_valid) ||
      any(abs(x$power - x$n_rejections / x$n_replicates) > 1e-10))
    stop("Invalid power or replicate counts in: ", file, call. = FALSE)

  data.frame(
    signal = ifelse(x$signal.method == "proportion", "Positive-part signal",
                    "Zero-part signal"),
    taxon = taxon,
    scenario = ifelse(x$scenario == "homogeneous", "homogeneous", "one_flipped"),
    signal_size = x$signal.size,
    method = x$method,
    source_method = x$source_method,
    n_replicates = x$n_replicates,
    n_valid = x$n_valid,
    n_rejected = x$n_rejections,
    failures = x$failures,
    calibration = x$calibration,
    n_boot = x$n_boot_requested,
    n_perm = x$n_perm_requested,
    source = basename(file),
    power = x$power,
    stringsAsFactors = FALSE
  )
}

plot_data <- do.call(rbind, lapply(taxa, read_power))
rownames(plot_data) <- NULL
if (any(plot_data$n_valid == 0L)) {
  warning("Some methods have no valid replicates; their power is NA.", call. = FALSE)
}

plot_data$signal <- factor(plot_data$signal,
                           levels = c("Positive-part signal", "Zero-part signal"))
plot_data$scenario <- factor(plot_data$scenario,
                             levels = c("homogeneous", "one_flipped"),
                             labels = c("Homogeneous scenario 1",
                                        "Heterogeneous scenario 2"))
plot_data$taxon <- factor(plot_data$taxon,
                          levels = c("Peptostreptococcaceae", "Enterobacteriaceae"))
plot_data$x <- unname(method_positions[plot_data$method])

# Match the existing colors and legend order; omit absent DM/ZIGDM families.
method_family <- function(method) {
  ifelse(grepl("MetaQCAT-O", method, fixed = TRUE), "MetaQCAT-O",
  ifelse(grepl("MetaQCAT*", method, fixed = TRUE), "MetaQCAT*",
  ifelse(grepl("MetaQCAT", method, fixed = TRUE), "MetaQCAT",
  ifelse(grepl("QCAT1", method, fixed = TRUE), "QCAT1", "QCAT2"))))
}
family_levels <- c("QCAT1", "QCAT2", "MetaQCAT", "MetaQCAT*", "MetaQCAT-O")
plot_data$family <- factor(method_family(plot_data$method), levels = family_levels)
palette <- c(
  "QCAT1" = "#009E8E",       # teal: QCAT.Cluster
  "QCAT2" = "#8B6BB1",       # purple: QCAT_GEE.Cluster
  "MetaQCAT" = "#E69F00",    # orange
  "MetaQCAT*" = "#D55E00",   # vermilion
  "MetaQCAT-O" = "#5C3317"   # deep brown
)

# Peptostreptococcaceae occupies the first row, Enterobacteriaceae the second.
# Within each row: positive homogeneous, zero homogeneous, positive flipped,
# zero flipped. Each plot gets its own natural y scale.
panel_specs <- expand.grid(
  taxon = levels(plot_data$taxon),
  scenario = levels(plot_data$scenario),
  signal = levels(plot_data$signal),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
panel_specs <- panel_specs[order(
  match(panel_specs$taxon, levels(plot_data$taxon)),
  match(panel_specs$scenario, levels(plot_data$scenario)),
  match(panel_specs$signal, levels(plot_data$signal))
), , drop = FALSE]

make_power_panel <- function(signal, scenario, taxon, show_y_title) {
  d <- plot_data[plot_data$signal == signal &
                   plot_data$scenario == scenario &
                   plot_data$taxon == taxon, , drop = FALSE]

  ggplot2::ggplot(d, ggplot2::aes(x = x, y = power, fill = family)) +
    ggplot2::geom_col(width = 0.78, colour = "black", linewidth = 0.42,
                      na.rm = TRUE) +
    ggplot2::scale_x_continuous(
      breaks = c(1.5, 4.5, 7.5, 11, 15),
      labels = c("P", "CD", "RD", "FE", "RE"),
      limits = c(0.35, 16.65),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = function(z) pretty(z, n = 4),
      labels = scales::label_number(accuracy = 0.1),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::scale_fill_manual(values = palette, breaks = family_levels,
                                drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      nrow = 1, override.aes = list(colour = "black", linewidth = 0.42)
    )) +
    ggplot2::labs(x = NULL, y = if (show_y_title) "Power" else NULL,
                  fill = NULL) +
    ggplot2::theme_classic(base_size = 17, base_family = "sans") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid = ggplot2::element_blank(),
      aspect.ratio = 1,
      panel.border = ggplot2::element_rect(fill = NA, colour = "black",
                                            linewidth = 0.8),
      axis.line = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", size = 22.95,
                                           colour = "#273444"),
      axis.text.y = ggplot2::element_text(face = "bold", size = 21.6,
                                           colour = "#273444"),
      axis.title.y = ggplot2::element_text(face = "bold", size = 19,
                                            colour = "#273444",
                                            margin = ggplot2::margin(r = 7)),
      legend.position = "bottom",
      legend.justification = "center",
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key.width = grid::unit(0.58, "cm"),
      legend.key.height = grid::unit(0.58, "cm"),
      legend.spacing.x = grid::unit(0.28, "cm"),
      legend.text = ggplot2::element_text(size = 15, colour = "#273444"),
      plot.margin = ggplot2::margin(10, 10, 8, 10)
    )
}

panels <- lapply(seq_len(nrow(panel_specs)), function(i) {
  make_power_panel(
    panel_specs$signal[i], panel_specs$scenario[i], panel_specs$taxon[i],
    show_y_title = i %in% c(1L, 5L)
  )
})
combined <- patchwork::wrap_plots(panels, ncol = 4, guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

# Add labels around the existing square panels, without changing their data,
# y scales, method order, black outlines, or collected bottom legend.
grDevices::pdf(file = NULL, width = 21, height = 12)
layout.device <- grDevices::dev.cur()
combined <- tryCatch(patchwork::patchworkGrob(combined),
  finally = grDevices::dev.off(layout.device))
panel.respect <- combined$respect
combined <- gtable::gtable_add_rows(
  combined, grid::unit(c(0.44, 0.38), "in"), pos = 0)
combined <- gtable::gtable_add_cols(
  combined, grid::unit(0.60, "in"), pos = 0)
if (is.matrix(panel.respect))
  combined$respect <- rbind(matrix(0, nrow = 2L, ncol = ncol(panel.respect) + 1L),
                           cbind(0, panel.respect))
panel.cells <- combined$layout[
  startsWith(combined$layout$name, "panel;") |
    grepl("^panel-[0-9]+$", combined$layout$name), , drop = FALSE]
panel.number <- as.integer(sub(".*-", "", panel.cells$name))
panel.cells <- panel.cells[match(seq_len(8L), panel.number), , drop = FALSE]
stopifnot(nrow(panel.cells) == 8L, !anyNA(panel.cells$l))

for (group in seq_len(2L)) {
  columns <- if (group == 1L) 1:2 else 3:4
  combined <- gtable::gtable_add_grob(combined,
    grid::textGrob(c("Low heterogeneity", "High heterogeneity")[group],
      gp = grid::gpar(fontfamily = "sans", fontface = "bold", fontsize = 28)),
    t = 1, b = 1, l = panel.cells$l[columns[1L]],
    r = panel.cells$r[columns[2L]], clip = "off",
    name = paste0("heterogeneity-header-", group))
}
for (column in seq_len(4L)) {
  combined <- gtable::gtable_add_grob(combined,
    grid::textGrob(rep(c("Small ZP change", "Large ZP change"), 2L)[column],
      gp = grid::gpar(fontfamily = "sans", fontface = "bold", fontsize = 25)),
    t = 2, b = 2, l = panel.cells$l[column], r = panel.cells$r[column],
    clip = "off", name = paste0("zp-header-", column))
}
for (row in seq_len(2L)) {
  cells <- (row - 1L) * 4L + seq_len(4L)
  combined <- gtable::gtable_add_grob(combined,
    grid::textGrob(taxa[row], rot = 90,
      gp = grid::gpar(fontfamily = "sans", fontface = "bold.italic",
                     fontsize = 27)),
    t = min(panel.cells$t[cells]), b = max(panel.cells$b[cells]),
    l = 1, r = 1, clip = "off", name = paste0("family-row-label-", row))
}

png_path <- file.path(plot_dir, "simu_res_cluster_power.png")
pdf_path <- file.path(plot_dir, "simu_res_cluster_power.pdf")
ggplot2::ggsave(png_path, combined, width = 21, height = 12.0,
                units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(pdf_path, combined, width = 21, height = 12.0,
                units = "in", device = grDevices::cairo_pdf, bg = "white")

csv_path <- file.path(plot_dir, "simu_res_cluster_power_plot_data.csv")
utils::write.csv(plot_data[c("signal", "taxon", "scenario", "method",
                             "source_method", "signal_size",
                             "n_replicates", "n_valid", "n_rejected", "power",
                             "failures", "calibration", "n_boot", "n_perm",
                             "source")],
                 csv_path, row.names = FALSE)
message("Saved 2 x 4 correlated power PNG: ", png_path)
message("Saved 2 x 4 correlated power PDF: ", pdf_path)
message("Saved plot data: ", csv_path)
