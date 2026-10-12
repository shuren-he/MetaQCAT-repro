# Plot MetaQCAT replicate p-values produced by CRC_*_Sim.R.
# Run from the project root: Rscript simulation/code/bar_plot.R
# Settings below select permutation or asymptotic results and the rejection cutoff.
# No simulation is run by this file. Only PDFs are saved to plot/simulation/.
# To load functions only: options(metaqcat.plot.source.only = TRUE); source(...)

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

input.dir <- file.path(project.root, "simulation", "outcome", "pos_part_res")
plot.dir <- file.path(project.root, "plot", "simulation")
inference <- "perm" # "perm" or "asy"
alpha <- 0.05

# Three MetaQCAT components in each FE/RE group.
method.positions <- c(
  "FE-MetaQCAT" = 1, "FE-MetaQCAT*" = 2, "FE-MetaQCAT-O" = 3,
  "RE-MetaQCAT" = 5, "RE-MetaQCAT*" = 6, "RE-MetaQCAT-O" = 7
)
method.colors <- c(
  MetaQCAT = "#E69F00", `MetaQCAT*` = "#D55E00", `MetaQCAT-O` = "#CC6677"
)
scenario.labels <- c(homogeneous = "Homogeneous",
                     two.null.studies = "Two null studies",
                     one.flipped.study = "One flipped study")

summarize_simulation_results <- function(input.dir, inference = "perm", alpha = 0.05) {
  inference <- match.arg(inference, c("perm", "asy"))
  stopifnot(length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)
  # Filename metadata matches file.tag in both current simulation scripts.
  pattern <- paste0(
    "^(Peptostreptococcaceae|Enterobacteriaceae)_(power|type1)_",
    "(proportion|count)_(homogeneous|two\\.null\\.studies|one\\.flipped\\.study)_",
    "([0-9.eE+-]+)_", inference, "_MetaQCAT_org2?\\.csv$"
  )
  files <- list.files(input.dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) return(NULL)
  do.call(rbind, lapply(files, function(file) {
    info <- regmatches(basename(file), regexec(pattern, basename(file)))[[1L]]
    values <- utils::read.csv(file, check.names = FALSE)
    methods <- sub(" \\((perm|asy)\\)$", "", names(values))
    # Also read older comparison CSVs, selecting only the MetaQCAT columns.
    keep <- grepl("^(FE|RE)-MetaQCAT(\\*|-O)?$", methods)
    values <- values[, keep, drop = FALSE]
    methods <- methods[keep]
    expected <- names(method.positions)
    if (inference == "asy") expected <- expected[!startsWith(expected, "RE-")]
    if (!nrow(values) || anyDuplicated(methods) || !setequal(methods, expected)) {
      stop("Unexpected method columns or no replicates in: ", file)
    }
    # read.csv() treats an entirely missing column as logical, not numeric.
    if (!all(vapply(values, function(x) is.numeric(x) || all(is.na(x)), logical(1)))) {
      stop("Non-numeric p-values in: ", file)
    }
    values <- as.matrix(data.frame(lapply(values, as.numeric)))
    if (any(!is.na(values) & (!is.finite(values) | values < 0 | values > 1))) {
      stop("P-values outside [0, 1] in: ", file)
    }
    n.valid <- colSums(!is.na(values))
    n.rejected <- colSums(values < alpha, na.rm = TRUE)
    rate <- rep(NA_real_, length(n.valid))
    rate[n.valid > 0] <- n.rejected[n.valid > 0] / n.valid[n.valid > 0]
    data.frame(
      lineage = info[2L], target = info[3L], signal = info[4L],
      scenario = info[5L], signal.size = as.numeric(info[6L]),
      inference = inference, alpha = alpha, method = methods,
      n.replicates = nrow(values), n.valid = as.integer(n.valid),
      n.failed = nrow(values) - as.integer(n.valid),
      n.rejected = as.integer(n.rejected), rejection.rate = rate,
      source.file = basename(file), row.names = NULL
    )
  }))
}

make_simulation_bar_plot <- function(d) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 to draw bar plots.")
  stopifnot(length(unique(d$target)) == 1L, length(unique(d$signal)) == 1L,
            length(unique(d$inference)) == 1L)
  d$x <- unname(method.positions[d$method])
  d$family <- factor(sub("^(FE|RE)-", "", d$method),
                     levels = names(method.colors))
  d$scenario <- factor(d$scenario, levels = names(scenario.labels),
                        labels = unname(scenario.labels))
  signal.label <- if (d$signal[1L] == "proportion") "Delta" else "p"
  panel <- paste0(d$lineage, "\n", signal.label, " = ", d$signal.size)
  panel.order <- order(match(d$lineage, c("Peptostreptococcaceae", "Enterobacteriaceae")),
                       d$signal.size)
  d$panel <- factor(panel, levels = unique(panel[panel.order]))
  y.label <- if (d$target[1L] == "power") "Power" else "Type I error rate"
  y.scale <- if (d$target[1L] == "power") {
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                               expand = ggplot2::expansion(mult = c(0, 0.04)))
  } else {
    ggplot2::scale_y_continuous(limits = c(0, NA),
                               expand = ggplot2::expansion(mult = c(0, 0.06)))
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = rejection.rate, fill = family)) +
    ggplot2::geom_col(width = 0.78, colour = "black", linewidth = 0.3, na.rm = TRUE) +
    ggplot2::facet_grid(rows = ggplot2::vars(scenario), cols = ggplot2::vars(panel)) +
    ggplot2::scale_x_continuous(
      breaks = if (d$inference[1L] == "perm") c(2, 6) else 2,
      labels = if (d$inference[1L] == "perm") c("FE", "RE") else "FE",
      expand = ggplot2::expansion(add = 0.8)
    ) +
    y.scale +
    ggplot2::scale_fill_manual(values = method.colors, drop = FALSE) +
    ggplot2::labs(x = NULL, y = y.label, fill = NULL,
      title = paste(if (d$target[1L] == "type1") "Null simulation" else
                      if (d$signal[1L] == "proportion") "Positive-part signal" else "Zero-part signal",
                    paste0("(", d$inference[1L], ")")),
      caption = paste0("Rejection cutoff = ", d$alpha[1L],
                       "; rates use valid replicates.")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(face = "bold"),
                   strip.text = ggplot2::element_text(face = "bold")) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))
  if (d$target[1L] == "type1") {
    p <- p + ggplot2::geom_hline(yintercept = d$alpha[1L], linetype = "dashed")
  }
  p
}

run_simulation_bar_plots <- function(input.dir, plot.dir, inference = "perm", alpha = 0.05) {
  summary <- summarize_simulation_results(input.dir, inference, alpha)
  if (is.null(summary)) {
    message("No matching simulation result CSVs in ", input.dir,
            ". Run CRC_Entero_Sim.R / CRC_Pepto_Sim.R first.")
    return(invisible(NULL))
  }
  if (any(summary$n.failed > 0L)) {
    warning("Some replicates failed; rates use only valid replicates.", call. = FALSE)
  }
  dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
  groups <- split(summary, interaction(summary$target, summary$signal, drop = TRUE))
  plots <- lapply(groups, function(d) {
    p <- make_simulation_bar_plot(d)
    tag <- paste(d$target[1L], d$signal[1L], inference, "bar_plot", sep = "_")
    width <- max(8, 6 * length(unique(paste(d$lineage, d$signal.size))))
    height <- 2.5 + 2.5 * length(unique(d$scenario))
    ggplot2::ggsave(file.path(plot.dir, paste0(tag, ".pdf")), p,
                    width = width, height = height, bg = "white")
    p
  })
  message("Saved simulation bar plot PDFs to ", plot.dir)
  invisible(list(summary = summary, plots = plots))
}

if (!isTRUE(getOption("metaqcat.plot.source.only", FALSE))) {
  run_simulation_bar_plots(input.dir, plot.dir, inference, alpha)
}
