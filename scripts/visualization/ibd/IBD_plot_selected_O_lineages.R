# Run from the MetaQCAT-repro repository root.
source(file.path("R", "project_paths.R"))

# Redraw selected IBD family summaries from the latest nonoverlapping analysis.
# First run IBD_extract_selected_O_lineages_latest.R to save MetaSub summaries.
# Use IBD_PLOT_STYLE=logratio for the descriptive count-ratio version.
# Use IBD_PLOT_LAYOUT=four_row to put all genera on one page while retaining
# both component summaries and their forest plots.
# Set IBD_SELECTED_DIR in both scripts to keep a new figure selection separate.

script.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
# Repository paths are defined by R/project_paths.R.

plot.style <- Sys.getenv("IBD_PLOT_STYLE", unset = "violin")
if (!plot.style %in% c("violin", "logratio"))
  stop("IBD_PLOT_STYLE must be 'violin' or 'logratio'.")
plot.layout <- Sys.getenv("IBD_PLOT_LAYOUT", unset = "full")
if (!plot.layout %in% c("full", "four_row", "two_row"))
  stop("IBD_PLOT_LAYOUT must be 'full', 'four_row', or 'two_row'.")
one.page.layout <- plot.layout %in% c("four_row", "two_row")
two.row.layout <- identical(plot.layout, "two_row")
required <- c("ggplot2", "forestplot", "gridExtra", "ggplotify")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

processed <- new.env(parent = emptyenv())
load(file.path(project.dir, "res", "IBD", "wild_score_bootstrap", "IBD_processed_bootstrap.RData"),
     envir = processed)

families <- c("f__Ruminococcaceae", "f__Peptostreptococcaceae")
requested.families <- commandArgs(trailingOnly = TRUE)
if (length(requested.families)) families <- requested.families
study.names <- names(processed$counts)
if (length(study.names) != 10L ||
    !identical(study.names, names(processed$design)) ||
    !identical(study.names, names(processed$subject_ids)))
  stop("Current IBD count, design and cluster lists are not aligned.")
study.labels <- setNames(paste0("IBD", seq_along(study.names)), study.names)
for (study in study.names) {
  if (!identical(rownames(processed$counts[[study]]),
                 rownames(processed$design[[study]])))
    stop("IBD sample rows are not aligned in ", study)
}

result.dir <- Sys.getenv("IBD_SELECTED_DIR", unset = "")
if (!nzchar(result.dir))
  result.dir <- file.path(project.dir, "outcome", "IBD", "selected_O_lineages")
if (!grepl("^/", result.dir)) result.dir <- file.path(project.dir, result.dir)
plot.dir <- file.path(result.dir, "plots")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  data.frame(plot_label = unname(study.labels),
             analysis_stratum = names(study.labels)),
  file.path(result.dir, "study_labels.csv"), row.names = FALSE
)

group.info <- do.call(rbind, lapply(study.names, function(study) {
  status <- processed$design[[study]][, "disease"]
  if (anyNA(status) || !all(status %in% c(0, 1)))
    stop("Disease indicator must be binary in ", study)
  data.frame(study = study.labels[[study]],
             Group = ifelse(status == 1, "IBD", "CTR"),
             stringsAsFactors = FALSE)
}))
group.info$study <- factor(group.info$study, levels = unname(study.labels))
group.info$Group <- factor(group.info$Group, levels = c("IBD", "CTR"))

plot.label <- function(taxon) {
  if (identical(taxon, "g__")) return("Unclassified")
  sub("^g__", "", taxon)
}

make_forest <- function(effects, taxon) {
  rows <- data.frame(study.id = unname(study.labels), mean = NA_real_,
                     lower = NA_real_, upper = NA_real_)
  observed <- effects[effects$taxon == taxon & effects$study != "Overall (FE)", ]
  if (nrow(observed)) {
    at <- match(observed$study, unname(study.labels))
    if (anyNA(at) || anyDuplicated(at)) stop("Unaligned study effects for ", taxon)
    rows$mean[at] <- observed$estimate
    rows$lower[at] <- observed$ci_lower
    rows$upper[at] <- observed$ci_upper
  }
  pooled <- effects[effects$taxon == taxon & effects$study == "Overall (FE)", ]
  if (nrow(pooled) != 1L) stop("Missing pooled FE effect for ", taxon)
  # forestplot uses the first row to establish its coordinate system; the
  # absent-study sentinel matches the CRC summary-plot convention.
  if (is.na(rows$mean[[1L]])) rows$mean[[1L]] <- -100
  forest <- forestplot::forestplot(
    rows, mean = mean, lower = lower, upper = upper, labeltext = study.id,
    # The seven-column version is substantially narrower than the original
    # multi-page layout. Scale only forest-plot typography to that cell size;
    # keep the estimates, intervals, symbols and axis limits unchanged.
    txt_gp = if (one.page.layout) {
      forestplot::fpTxtGp(
        label = grid::gpar(cex = 0.48),
        summary = grid::gpar(cex = 0.52),
        ticks = grid::gpar(cex = 0.48),
        xlab = grid::gpar(cex = 0.52)
      )
    } else {
      forestplot::fpTxtGp(xlab = grid::gpar(cex = 0.75))
    },
    xlab = "one-step score summaries", ci.vertices = TRUE,
    ci.vertices.height = 0.1
  ) |>
    forestplot::fp_set_style(box = "black", line = "black",
                             summary = "royalblue", align = "l",
                             hrz_lines = "#999999") |>
    forestplot::fp_add_lines() |>
    forestplot::fp_add_header(study.id = "Study") |>
    forestplot::fp_append_row(
      mean = pooled$estimate, lower = pooled$ci_lower,
      upper = pooled$ci_upper, study.id = "Overall (FE)", is.summary = TRUE
    )
  ggplotify::grid2grob(print(forest))
}

empty.panel <- function() ggplot2::ggplot() + ggplot2::theme_void()

extract.legend <- function(plot) {
  built <- ggplot2::ggplotGrob(plot)
  at <- grep("^guide-box", built$layout$name)
  at <- at[!vapply(built$grobs[at], inherits, logical(1), "zeroGrob")]
  if (length(at) != 1L) stop("Expected one visible IBD/CTR legend")
  built$grobs[[at]]
}

plot.family <- function(family, summary, file.stem) {
  child.counts <- summary$counts
  rarefied <- summary$rarefied_counts
  if (!identical(names(child.counts), study.names) ||
      !identical(names(rarefied), study.names)) {
    stop("MetaSub study order differs from processed IBD data for ", family)
  }
  for (study in study.names) {
    if (nrow(child.counts[[study]]) != nrow(processed$counts[[study]]) ||
        nrow(rarefied[[study]]) != nrow(processed$counts[[study]]) ||
        !identical(rownames(child.counts[[study]]),
                   rownames(processed$counts[[study]])) ||
        !identical(rownames(rarefied[[study]]),
                   rownames(processed$counts[[study]]))) {
      stop("MetaSub sample rows differ in ", family, " / ", study)
    }
  }
  Y <- do.call(rbind, child.counts)
  Y.zero <- do.call(rbind, rarefied)
  taxa <- colnames(Y)
  if (!setequal(taxa, colnames(Y.zero)) || anyDuplicated(taxa))
    stop("MetaSub positive/zero taxa differ for ", family)
  Y.zero <- Y.zero[, taxa, drop = FALSE]
  reference <- summary$positive$reference.taxon
  if (!reference %in% taxa)
    stop("MetaSub reference taxon is missing for ", family)
  taxon.order <- taxa[order(colSums(Y), decreasing = TRUE)]
  taxon.order <- c(taxon.order[taxon.order != reference], reference)
  if (anyDuplicated(vapply(taxon.order, plot.label, character(1))))
    stop("Plot genus labels are not unique for ", family)

  positive.effects <- summary$positive$estimates
  zero.effects <- summary$zero$estimates
  remap.studies <- function(effects) {
    real <- effects$study != "Overall (FE)"
    if (!all(effects$study[real] %in% study.names))
      stop("Unknown IBD stratum in MetaSub effects for ", family)
    effects$study[real] <- unname(study.labels[effects$study[real]])
    effects
  }
  positive.effects <- remap.studies(positive.effects)
  zero.effects <- remap.studies(zero.effects)
  required.columns <- c("study", "taxon", "estimate", "se", "ci_lower", "ci_upper")
  if (!all(required.columns %in% names(positive.effects)) ||
      !all(required.columns %in% names(zero.effects))) {
    stop("MetaSub effect table columns are incomplete for ", family)
  }
  positive.selected <- unique(positive.effects$taxon[
    positive.effects$study == "Overall (FE)"])
  zero.selected <- unique(zero.effects$taxon[
    zero.effects$study == "Overall (FE)"])

  # The log-ratio layer is descriptive. The forest below remains the score
  # estimate, which uses all family counts and is not a mean of these dots.
  keep <- rowSums(Y) > 0
  position <- which(keep)
  positive.data <- do.call(rbind, lapply(taxon.order, function(taxon) {
    numerator <- Y[position, taxon]
    ref.count <- Y[position, reference]
    both.zero <- numerator == 0 & ref.count == 0
    value <- if (plot.style == "logratio") {
      log((numerator + 0.5) / (ref.count + 0.5))
    } else {
      numerator / (numerator + ref.count)
    }
    value[both.zero] <- NA_real_
    value[!is.finite(value)] <- NA_real_
    data.frame(taxon = taxon, study = group.info$study[position],
               Group = group.info$Group[position], value = value,
               numerator = numerator, ref.count = ref.count,
               family.depth = rowSums(Y)[position], both.zero = both.zero)
  }))
  estimated <- positive.effects[positive.effects$study != "Overall (FE)",
                                c("study", "taxon")]
  positive.key <- paste(estimated$study, estimated$taxon, sep = "\r")
  positive.data$model.valid <- paste(positive.data$study, positive.data$taxon,
                                     sep = "\r") %in% positive.key
  positive.data$value[!positive.data$model.valid] <- NA_real_

  zero.data <- do.call(rbind, lapply(taxon.order, function(taxon) {
    data.frame(taxon = taxon, study = group.info$study,
               Group = group.info$Group,
               value = as.numeric(Y.zero[, taxon] == 0))
  }))
  estimated <- zero.effects[zero.effects$study != "Overall (FE)",
                            c("study", "taxon")]
  zero.key <- paste(estimated$study, estimated$taxon, sep = "\r")
  zero.data$value[!(paste(zero.data$study, zero.data$taxon,
                          sep = "\r") %in% zero.key)] <- NA_real_
  zero.summary <- aggregate(value ~ taxon + study + Group,
                            zero.data[is.finite(zero.data$value), ], mean)

  has.positive <- vapply(taxon.order, function(taxon)
    taxon %in% positive.selected &&
      any(is.finite(positive.data$value[positive.data$taxon == taxon])),
    logical(1))
  has.zero <- vapply(taxon.order, function(taxon)
    taxon %in% zero.selected &&
      any(zero.summary$taxon == taxon), logical(1))
  # Spread a large family evenly across pages of at most five taxa; in a
  # one-page layout, place all taxa on one page.
  # Keep one legend at the last nonempty chart of each component on each page.
  n.taxa <- length(taxon.order)
  n.pages <- if (one.page.layout) 1L else ceiling(n.taxa / 5L)
  page.sizes <- if (one.page.layout) n.taxa else rep(n.taxa %/% n.pages, n.pages)
  if (!one.page.layout && n.taxa %% n.pages)
    page.sizes[seq_len(n.taxa %% n.pages)] <-
      page.sizes[seq_len(n.taxa %% n.pages)] + 1L
  page.indices <- split(seq_len(n.taxa), rep(seq_len(n.pages), page.sizes))
  page.legend <- function(plottable) {
    answer <- vapply(page.indices, function(index) {
      available <- taxon.order[index[plottable[index]]]
      if (length(available)) tail(available, 1L) else NA_character_
    }, character(1))
    unname(answer[!is.na(answer)])
  }
  positive.legend <- page.legend(has.positive)
  zero.legend <- page.legend(has.zero)
  if (!length(positive.legend) || !length(zero.legend))
    stop("No positive or zero component is plottable for ", family)

  palette <- c(IBD = "#CD3333", CTR = "steelblue")
  quartiles <- function(x) c(ymin = stats::quantile(x, .25, na.rm = TRUE),
                             ymax = stats::quantile(x, .75, na.rm = TRUE))
  # Keep each chart in the same-sized cell. A legend gets its own narrow
  # sidecar beside the last nonempty chart, never reducing the chart width.
  chart.width <- if (two.row.layout) 1.5 else if (one.page.layout) 1.35 else 4.0
  legend.width <- if (two.row.layout) 0.55 else if (one.page.layout) 0.62 else 0.9
  chart.inset <- if (two.row.layout) 0.03 else if (one.page.layout) 0.05 else 0.1
  column.gap <- if (one.page.layout) 0.08 else 0.22
  axis.title.width <- if (two.row.layout) 0.4 else if (one.page.layout) 0.55 else 0.62
  heading.height <- if (two.row.layout) 0.28 else if (one.page.layout) 0.38 else 0.55
  chart.height <- if (two.row.layout) 1.55 else if (one.page.layout) 1.45 else 3.75
  forest.height <- if (two.row.layout) 0 else if (one.page.layout) 2.15 else 4.1
  # ggplot mark sizes are physical units and do not shrink when a full-size
  # chart is fitted into the seven-column layout. Scale marks with the panel.
  mark.scale <- if (one.page.layout) chart.width / 4 else 1
  row.heights <- if (two.row.layout) {
    grid::unit(c(heading.height, chart.height, chart.height), "in")
  } else {
    grid::unit(c(heading.height, chart.height, forest.height,
                 chart.height, forest.height), "in")
  }
  panels <- vector("list", length(taxon.order))
  column.widths <- rep(chart.width, length(taxon.order))
  for (j in seq_along(taxon.order)) {
    taxon <- taxon.order[[j]]
    pos.legend.grob <- grid::nullGrob()
    zero.legend.grob <- grid::nullGrob()
    if (has.positive[[j]]) {
      d <- positive.data[positive.data$taxon == taxon, ]
      if (plot.style == "logratio") {
        pooled <- stats::aggregate(
          cbind(numerator, ref.count) ~ study + Group,
          data = d[d$model.valid, , drop = FALSE], FUN = sum)
        pooled$value <- ifelse(pooled$numerator > 0 & pooled$ref.count > 0,
                               log(pooled$numerator / pooled$ref.count),
                               NA_real_)
        pooled <- pooled[is.finite(pooled$value), , drop = FALSE]
        pooled$x <- as.numeric(pooled$study) +
          ifelse(pooled$Group == "IBD", -.18, .18)
        crc <- pooled[pooled$Group == "IBD", c("study", "value")]
        ctr <- pooled[pooled$Group == "CTR", c("study", "value")]
        paired <- merge(crc, ctr, by = "study", suffixes = c(".IBD", ".CTR"))
        links <- data.frame(x = as.numeric(paired$study) - .18,
                            xend = as.numeric(paired$study) + .18,
                            y = paired$value.IBD, yend = paired$value.CTR)
        dots <- d[is.finite(d$value), , drop = FALSE]
        dots$x <- as.numeric(dots$study) +
          ifelse(dots$Group == "IBD", -.18, .18)
        # Cap only marker size, not the plotted count ratio or pooled totals.
        dots$depth.size <- pmin(log1p(dots$family.depth),
                                 stats::quantile(log1p(dots$family.depth),
                                                 .95, na.rm = TRUE))
        violin <- ggplot2::ggplot() +
          ggplot2::geom_point(
            data = dots,
            ggplot2::aes(x = x, y = value, colour = Group, size = depth.size),
            position = ggplot2::position_jitter(width = .06, height = 0,
                                                seed = 1234L),
            alpha = .24, show.legend = FALSE) +
          ggplot2::geom_segment(
            data = links,
            ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            linewidth = .65 * mark.scale, colour = "#444444") +
          ggplot2::geom_point(
            data = pooled,
            ggplot2::aes(x = x, y = value, fill = Group),
            shape = 21, size = 3.2 * mark.scale,
            stroke = .7 * mark.scale, colour = "black") +
          ggplot2::scale_x_continuous(
            breaks = seq_along(study.names), labels = unname(study.labels),
            limits = c(.5, length(study.names) + .5)) +
          ggplot2::scale_size_continuous(
            range = c(.45, 2.2) * mark.scale, guide = "none") +
          ggplot2::scale_colour_manual(values = palette, guide = "none") +
          ggplot2::scale_fill_manual(values = palette) +
          ggplot2::theme_bw(base_size = if (two.row.layout) 8 else 11) +
          ggplot2::labs(y = NULL, x = NULL) +
          ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.key.size = grid::unit(
              if (one.page.layout) .4 else .75, "lines"),
            legend.text = ggplot2::element_text(
              size = if (one.page.layout) 6 else 9),
            aspect.ratio = 1,
            axis.text.x = ggplot2::element_text(
              angle = 90, hjust = 1,
              size = if (two.row.layout) 6.5 else
                if (one.page.layout) 5.5 else ggplot2::rel(1)),
            axis.text.y = ggplot2::element_text(
              size = if (two.row.layout) 6.5 else
                if (one.page.layout) 6 else ggplot2::rel(1)),
            plot.margin = ggplot2::margin(
              .1, 0, .1,
              if (two.row.layout) .15 else
                if (one.page.layout) .2 else .6, "cm")
          )
      } else {
        d <- d[is.finite(d$value), , drop = FALSE]
        violin <- ggplot2::ggplot(d, ggplot2::aes(x = study, y = value,
                                                  fill = Group)) +
          ggplot2::geom_violin(trim = TRUE, scale = "width",
                               quantile.linetype = 1, quantiles = 0.5,
                               width = .8,
                               position = ggplot2::position_dodge(width = .9),
                               na.rm = TRUE) +
          ggplot2::stat_summary(fun = quartiles, geom = "line",
                                linewidth = .5 * mark.scale,
                                position = ggplot2::position_dodge(width = .9),
                                color = "white", alpha = .5,
                                show.legend = FALSE, na.rm = TRUE) +
          ggplot2::stat_summary(fun = "mean", geom = "point", shape = 2,
                                size = 1.5 * mark.scale,
                                position = ggplot2::position_dodge(width = .9),
                                show.legend = FALSE, na.rm = TRUE) +
          ggplot2::scale_x_discrete(drop = FALSE) +
          ggplot2::scale_fill_manual(values = palette) +
          ggplot2::theme_bw(base_size = if (two.row.layout) 8 else 11) +
          ggplot2::labs(y = NULL, x = NULL) +
          ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.key.size = grid::unit(
              if (one.page.layout) .4 else .75, "lines"),
            legend.text = ggplot2::element_text(
              size = if (one.page.layout) 6 else 9),
            aspect.ratio = 1,
            axis.text.x = ggplot2::element_text(
              angle = 90, hjust = 1,
              size = if (two.row.layout) 6.5 else
                if (one.page.layout) 5.5 else ggplot2::rel(1)),
            axis.text.y = ggplot2::element_text(
              size = if (two.row.layout) 6.5 else
                if (one.page.layout) 6 else ggplot2::rel(1)),
            plot.margin = ggplot2::margin(
              .1, 0, .1,
              if (two.row.layout) .15 else
                if (one.page.layout) .2 else .6, "cm")
          )
      }
      if (taxon %in% positive.legend)
        pos.legend.grob <- extract.legend(violin)
      violin <- violin + ggplot2::theme(legend.position = "none")
      violin <- gridExtra::arrangeGrob(
        grid::nullGrob(), violin, nrow = 1,
        widths = grid::unit(c(chart.inset, chart.width - chart.inset), "in"))
      pos.forest <- if (two.row.layout) empty.panel() else
        make_forest(positive.effects, taxon)
    } else {
      violin <- empty.panel()
      pos.forest <- empty.panel()
    }
    if (has.zero[[j]]) {
      d <- zero.summary[zero.summary$taxon == taxon, ]
      bars <- ggplot2::ggplot(d, ggplot2::aes(x = study, y = value,
                                              fill = Group)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = .8),
                          width = .7, colour = "black", na.rm = TRUE) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::scale_fill_manual(values = palette) +
        ggplot2::theme_bw(base_size = if (two.row.layout) 8 else 11) +
        ggplot2::labs(y = NULL, x = NULL) +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.key.size = grid::unit(
            if (one.page.layout) .4 else .75, "lines"),
          legend.text = ggplot2::element_text(
            size = if (one.page.layout) 6 else 9),
          aspect.ratio = 1,
          axis.text.x = ggplot2::element_text(
            angle = 90, hjust = 1,
            size = if (two.row.layout) 6.5 else
              if (one.page.layout) 5.5 else ggplot2::rel(1)),
          axis.text.y = ggplot2::element_text(
            size = if (two.row.layout) 6.5 else
              if (one.page.layout) 6 else ggplot2::rel(1)),
          plot.margin = ggplot2::margin(.1, 0, .1,
                                        if (two.row.layout) .05 else .2, "cm")
        )
      if (taxon %in% zero.legend)
        zero.legend.grob <- extract.legend(bars)
      bars <- bars + ggplot2::theme(legend.position = "none")
      bars <- gridExtra::arrangeGrob(
        grid::nullGrob(), bars, nrow = 1,
        widths = grid::unit(c(chart.inset, chart.width - chart.inset), "in"))
      zero.forest <- if (two.row.layout) empty.panel() else
        make_forest(zero.effects, taxon)
    } else {
      bars <- empty.panel()
      zero.forest <- empty.panel()
    }
    title <- paste0(plot.label(taxon),
                    if (identical(taxon, reference)) " (Ref.)" else "")
    title.grob <- grid::textGrob(
      title, x = grid::unit(.55, "npc"),
      gp = grid::gpar(fontsize = if (two.row.layout) 10 else
                        if (one.page.layout) 11 else 18,
                      fontface = "bold"))
    chart.grobs <- if (two.row.layout) {
      list(title.grob, violin, bars)
    } else {
      list(title.grob, violin, pos.forest, bars, zero.forest)
    }
    chart.column <- gridExtra::arrangeGrob(
      grobs = chart.grobs, ncol = 1, heights = row.heights)
    if (taxon %in% positive.legend || taxon %in% zero.legend) {
      legend.column <- if (two.row.layout) {
        gridExtra::arrangeGrob(grid::nullGrob(), pos.legend.grob,
                               zero.legend.grob, ncol = 1,
                               heights = row.heights)
      } else {
        gridExtra::arrangeGrob(
          grid::nullGrob(), pos.legend.grob, grid::nullGrob(),
          zero.legend.grob, grid::nullGrob(),
          ncol = 1, heights = row.heights
        )
      }
      panels[[j]] <- gridExtra::arrangeGrob(
        chart.column, legend.column, nrow = 1,
        widths = grid::unit(c(chart.width, legend.width), "in")
      )
      column.widths[[j]] <- chart.width + legend.width
    } else {
      panels[[j]] <- chart.column
    }
  }

  positive.axis.title <- if (plot.style == "logratio") {
    "Log count ratio (taxon/ref.)"
  } else {
    "Taxon proportion relative to Ref."
  }
  if (two.row.layout) positive.axis.title <- "Positive part"
  zero.axis.title <- if (two.row.layout) "Zero proportion" else
    "Proportion of zero counts"
  pages <- lapply(page.indices, function(index) {
    axis.titles <- gridExtra::arrangeGrob(
      grid::nullGrob(),
      grid::textGrob(positive.axis.title, rot = 90,
                     gp = grid::gpar(fontsize = if (two.row.layout) 9 else
                                       if (one.page.layout) 8 else 16,
                                     fontface = "bold")),
      grid::textGrob(zero.axis.title, rot = 90,
                     gp = grid::gpar(fontsize = if (two.row.layout) 9 else
                                       if (one.page.layout) 8 else 16,
                                     fontface = "bold")),
      ncol = 1,
      heights = if (two.row.layout) {
        grid::unit(c(heading.height, chart.height, chart.height), "in")
      } else {
        grid::unit(c(heading.height, chart.height + forest.height,
                     chart.height + forest.height), "in")
      })
    page.grobs <- list(axis.titles)
    page.columns <- axis.title.width
    for (k in seq_along(index)) {
      page.grobs[[length(page.grobs) + 1L]] <- panels[[index[[k]]]]
      page.columns <- c(page.columns, column.widths[[index[[k]]]])
      if (k < length(index)) {
        page.grobs[[length(page.grobs) + 1L]] <- grid::nullGrob()
        page.columns <- c(page.columns, column.gap)
      }
    }
    gridExtra::arrangeGrob(grobs = page.grobs, nrow = 1,
                           widths = grid::unit(page.columns, "in"))
  })
  page.widths <- vapply(page.indices, function(index)
    axis.title.width + sum(column.widths[index]) +
      column.gap * (length(index) - 1L), numeric(1))
  page.margin <- 0.2
  figure.width <- max(page.widths)
  figure.height <- if (two.row.layout) {
    sum(c(heading.height, chart.height, chart.height))
  } else {
    sum(c(heading.height, chart.height, forest.height,
          chart.height, forest.height))
  }
  draw.page <- function(index, y = grid::unit(.5, "npc")) {
    grid::pushViewport(grid::viewport(y = y,
                                      width = grid::unit(page.widths[[index]], "in"),
                                      height = grid::unit(figure.height, "in")))
    grid::grid.draw(pages[[index]])
    grid::upViewport()
  }
  pdf.path <- file.path(plot.dir, paste0(file.stem, ".pdf"))
  grDevices::pdf(pdf.path, width = figure.width + 2 * page.margin,
                 height = figure.height + 2 * page.margin)
  for (index in seq_along(pages)) {
    grid::grid.newpage()
    draw.page(index)
  }
  grDevices::dev.off()
  png.path <- file.path(plot.dir, paste0(file.stem, ".png"))
  grDevices::png(png.path, width = figure.width + 2 * page.margin,
                 height = length(pages) * figure.height + 2 * page.margin,
                 units = "in", res = 150, type = "cairo")
  grid::grid.newpage()
  for (index in seq_along(pages)) {
    y <- grid::unit(page.margin +
                     (length(pages) - index + .5) * figure.height, "in")
    draw.page(index, y = y)
  }
  grDevices::dev.off()
  message("Saved ", pdf.path)
  invisible(c(pdf = pdf.path, png = png.path))
}

for (family in families) {
  stem <- sub("^f__", "", family)
  summary.path <- file.path(result.dir, paste0(stem, "_latest_summary.rds"))
  if (!file.exists(summary.path))
    stop("Run IBD_extract_selected_O_lineages_latest.R first: ", summary.path)
  summary <- readRDS(summary.path)
  if (!identical(summary$lineage, family) ||
      !identical(summary$min.prev, 0.1) ||
      !identical(summary$min.depth, 0) ||
      !identical(summary$seed, 1234L))
    stop("Saved summary does not match current IBD settings for ", family)
  output.stem <- paste0(stem,
                        if (plot.style == "logratio") "_logratio" else "",
                        if (two.row.layout) "_two_row" else
                          if (one.page.layout) "_four_row" else "")
  plot.family(family, summary, output.stem)
}
