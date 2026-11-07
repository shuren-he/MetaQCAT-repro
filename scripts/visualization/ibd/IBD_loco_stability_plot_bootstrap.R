# Run from the MetaQCAT-repro repository root.
source(file.path("R", "project_paths.R"))

# IBD leave-one-cohort-out heatmap for the current bootstrap results.
# Run from the repository root: Rscript IBD_loco_stability_plot_bootstrap.R
library(ggplot2)
library(patchwork)

script.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
# Repository paths are defined by R/project_paths.R.
result.file <- file.path(script.dir, "res", "IBD", "wild_score_bootstrap",
                         "leave_one_cohort_out", "IBD_loso_lineage_results.csv")
plot.dir <- file.path(script.dir, "outcome", "IBD")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
results <- read.csv(result.file, stringsAsFactors = FALSE, check.names = FALSE)

methods <- c("FE-MetaQCAT-O", "RE-MetaQCAT-O")
cohorts <- c("CS-PRISM", "Herfarth_CCFA_Microbiome_3B_combined", "HMP2",
             "Jansson_Lamendella_Crohns", "MucosalIBD", "Pouchitis", "RISK")
all.fits <- c("Full_data_300", "Full_data_1000", cohorts)
stopifnot(all(c("method", "omitted_cohort", "lineage", "q") %in% names(results)),
          !anyDuplicated(results[c("method", "omitted_cohort", "lineage")]),
          setequal(results$method, methods),
          setequal(results$omitted_cohort, all.fits),
          all(table(results$method, results$omitted_cohort) == 45L))

# Match the full-data baseline to the 300-draw cohort-omission fits.
# The primary 1,000-draw analysis is reported separately, not in this plot.
results <- results[results$omitted_cohort != "Full_data_1000", , drop = FALSE]
display.fits <- c("Full_data_300", cohorts)
display.labels <- c("All seven cohorts", "Drop CS-PRISM", "Drop Herfarth",
                    "Drop HMP2", "Drop Jansson--Lamendella", "Drop MucosalIBD",
                    "Drop Pouchitis", "Drop RISK")

panels <- list()
plot.data <- list()
for (method in methods) {
  d <- results[results$method == method, , drop = FALSE]
  discoveries <- tapply(d$q < 0.05, d$lineage, sum, na.rm = TRUE)
  discoveries <- discoveries[discoveries > 0]
  lineage.order <- names(discoveries)[order(discoveries, names(discoveries))]
  d <- d[d$lineage %in% lineage.order, , drop = FALSE]
  d$lineage <- factor(d$lineage, levels = rev(lineage.order))
  d$omitted_cohort <- factor(d$omitted_cohort, levels = display.fits,
                             labels = display.labels)
  d$status <- factor(ifelse(is.na(d$q), "Not estimable",
    ifelse(d$q < 0.05, "q < 0.05", "Not significant")),
    levels = c("q < 0.05", "Not significant", "Not estimable"))
  plot.data[[method]] <- d

  panels[[method]] <- ggplot(d, aes(omitted_cohort, lineage, fill = status)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_manual(values = c("q < 0.05" = "#00858E",
      "Not significant" = "#E4E9ED", "Not estimable" = "#F1E8D5"),
      drop = FALSE) +
    scale_x_discrete(expand = c(0, 0), drop = FALSE) +
    scale_y_discrete(labels = function(x) sub("^f__", "", x),
                     expand = c(0, 0)) +
    labs(x = NULL, y = NULL, fill = NULL,
         title = if (method == methods[1L]) "FE" else "RE") +
    theme_minimal(base_size = 12, base_family = "sans") +
    theme(panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, face = "bold",
                                 colour = "#4A4A4A"),
      axis.text.y = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = "bottom", legend.text = element_text(size = 12),
      legend.key.width = grid::unit(0.55, "cm"),
      legend.key.height = grid::unit(0.55, "cm"),
      plot.margin = margin(8, 8, 8, 8))
}

combined <- wrap_plots(panels, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")
pdf.file <- file.path(plot.dir, "IBD_loco_stability_bootstrap.pdf")
png.file <- file.path(plot.dir, "IBD_loco_stability_bootstrap.png")
ggsave(pdf.file, combined, width = 13, height = 637 / 72,
       units = "in", device = grDevices::cairo_pdf, bg = "white")
ggsave(png.file, combined, width = 13, height = 637 / 72,
       units = "in", dpi = 300, bg = "white")
write.csv(do.call(rbind, plot.data), file.path(plot.dir,
  "IBD_loco_stability_bootstrap_plot_data.csv"), row.names = FALSE)
message("Saved: ", pdf.file)
message("Saved: ", png.file)
