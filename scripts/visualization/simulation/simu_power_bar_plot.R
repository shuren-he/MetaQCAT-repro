# Simulation power: all methods evaluated on shared generated datasets.
# Run from the repository root:
# Rscript scripts/visualization/simulation/simu_power_bar_plot.R

required <- c("ggplot2", "patchwork", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

source("R/project_paths.R", local = TRUE)
project_dir <- project.dir

input_dir <- file.path(project_dir, "res")
plot_dir <- file.path(project_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Preserve the paper's method order and group spacing.
method_positions <- c(
  "P-DM (perm)" = 1,
  "P-ZIGDM (perm)" = 2,
  "P-QCAT1 (perm)" = 3,
  "P-QCAT2 (perm)" = 4,
  "CD-DM (perm)" = 6,
  "CD-ZIGDM (perm)" = 7,
  "CD-QCAT1 (perm)" = 8,
  "CD-QCAT2 (perm)" = 9,
  "RD-DM (perm)" = 11,
  "RD-ZIGDM (perm)" = 12,
  "RD-QCAT1 (perm)" = 13,
  "RD-QCAT2 (perm)" = 14,
  "FE-MetaQCAT (boot)" = 16,
  "FE-MetaQCAT* (boot)" = 17,
  "FE-MetaQCAT-O (boot)" = 18,
  "RE-MetaQCAT (boot)" = 20,
  "RE-MetaQCAT* (boot)" = 21,
  "RE-MetaQCAT-O (boot)" = 22
)

taxa <- c("Peptostreptococcaceae", "Enterobacteriaceae")
read_power <- function(taxon) {
  file <- file.path(input_dir, paste0(taxon, "_power_summary.csv"))
  if (!file.exists(file)) stop("Missing simulation summary: ", file, call. = FALSE)
  x <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("lineage", "simulation_target", "scenario", "signal.method",
              "method", "rejection_rate", "n_rejections", "n_replicates",
              "n_valid", "failures", "n_perm_configured", "calibration",
              "n_boot_requested")
  if (!all(needed %in% names(x)) || anyDuplicated(names(x)))
    stop("Unexpected summary columns in: ", file, call. = FALSE)
  if (any(x$lineage != taxon) || any(x$simulation_target != "power"))
    stop("Unexpected lineage or simulation target in: ", file, call. = FALSE)

  is.meta <- grepl("MetaQCAT", x$method, fixed = TRUE)
  if (any(x$calibration[is.meta] != "wild-score bootstrap") ||
      any(x$calibration[!is.meta] != "permutation"))
    stop("Unexpected method calibration in: ", file, call. = FALSE)
  x$method[is.meta] <- paste0(x$method[is.meta], " (boot)")
  x$source <- basename(file)

  expected <- expand.grid(
    method = names(method_positions),
    signal.method = c("proportion", "count"),
    scenario = c("homogeneous", "one.flipped.study"),
    stringsAsFactors = FALSE
  )
  key <- function(d) paste(d$method, d$signal.method, d$scenario, sep = "\r")
  if (nrow(x) != nrow(expected) || anyDuplicated(key(x)) ||
      !setequal(key(x), key(expected)))
    stop("Expected 18 methods in each of four settings in: ", file,
         call. = FALSE)
  x <- x[match(key(expected), key(x)), , drop = FALSE]
  # Every method uses the planned replicate denominator; failures are retained.
  denominator <- x$n_replicates
  if (any(!is.finite(x$rejection_rate) | x$rejection_rate < 0 |
          x$rejection_rate > 1 | denominator < 1L | x$n_valid < 0L |
          x$n_valid > x$n_replicates | x$n_rejections < 0L |
          x$n_rejections > x$n_valid |
          x$failures != x$n_replicates - x$n_valid) ||
      any(abs(x$rejection_rate - x$n_rejections / denominator) > 1e-10))
    stop("Invalid resampling power or replicate counts in: ", file,
         call. = FALSE)

  data.frame(
    signal = ifelse(x$signal.method == "proportion", "Positive-part signal",
                    "Zero-part signal"),
    taxon = taxon,
    scenario = ifelse(x$scenario == "homogeneous", "homogeneous", "one_flipped"),
    method = x$method,
    n_replicates = x$n_replicates,
    n_valid = x$n_valid,
    n_rejected = x$n_rejections,
    failures = x$failures,
    calibration = x$calibration,
    n_boot = x$n_boot_requested,
    n_perm = x$n_perm_configured,
    source = x$source,
    power = x$rejection_rate,
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

# Method identity is carried by fill color throughout the figure. Three warm
# hues denote MetaQCAT components; four cool hues denote comparators.
method_family <- function(method) {
  ifelse(grepl("MetaQCAT-O", method, fixed = TRUE), "MetaQCAT-O",
  ifelse(grepl("MetaQCAT*", method, fixed = TRUE), "MetaQCAT*",
  ifelse(grepl("MetaQCAT", method, fixed = TRUE), "MetaQCAT",
  ifelse(grepl("QCAT1", method, fixed = TRUE), "QCAT1",
  ifelse(grepl("QCAT2", method, fixed = TRUE), "QCAT2",
  ifelse(grepl("ZIGDM", method, fixed = TRUE), "ZIGDM", "DM"))))))
}
family_levels <- c("DM", "ZIGDM", "QCAT1", "QCAT2",
                   "MetaQCAT", "MetaQCAT*", "MetaQCAT-O")
plot_data$family <- factor(method_family(plot_data$method), levels = family_levels)
palette <- c(
  "DM" = "#3B5BA9",          # deep blue
  "ZIGDM" = "#56B4E9",       # sky blue
  "QCAT1" = "#009E8E",       # teal
  "QCAT2" = "#8B6BB1",       # purple
  "MetaQCAT" = "#E69F00",    # orange
  "MetaQCAT*" = "#D55E00",   # vermilion
  "MetaQCAT-O" = "#CC6677"   # rose
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
      breaks = c(2.5, 7.5, 12.5, 17, 21),
      labels = c("P", "CD", "RD", "FE", "RE"),
      limits = c(0.35, 22.65),
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
      panel.border = ggplot2::element_rect(fill = NA, colour = "black",
                                            linewidth = 0.8),
      axis.line = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", size = 17,
                                           colour = "#273444"),
      axis.text.y = ggplot2::element_text(face = "bold", size = 16,
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

png_path <- file.path(plot_dir, "simu_res_power.png")
pdf_path <- file.path(plot_dir, "simu_res_power.pdf")
ggplot2::ggsave(png_path, combined, width = 21, height = 11.1,
                units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(pdf_path, combined, width = 21, height = 11.1,
                units = "in", device = grDevices::cairo_pdf, bg = "white")

csv_path <- file.path(plot_dir, "simu_res_power_plot_data.csv")
utils::write.csv(plot_data[c("signal", "taxon", "scenario", "method",
                             "n_replicates", "n_valid", "n_rejected", "power",
                             "failures", "calibration", "n_boot", "n_perm",
                             "source")],
                 csv_path, row.names = FALSE)
message("Saved 2 x 4 power PNG: ", png_path)
message("Saved 2 x 4 power PDF: ", pdf_path)
message("Saved plot data: ", csv_path)
