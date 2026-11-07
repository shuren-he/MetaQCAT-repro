# Run from the MetaQCAT-repro repository root.
source(file.path("R", "project_paths.R"))

# Recreate the CRC four-row summary plots. Set CRC_PLOT_STYLE=logratio for the
# descriptive log-ratio version, and CRC_REUSE_SUMMARY=true to redraw the
# violin version without refitting the stored MetaQCAT summary.
# Example: CRC_REUSE_SUMMARY=true Rscript CRC_plot_selected_O_lineages.R
#          CRC_PLOT_STYLE=logratio Rscript CRC_plot_selected_O_lineages.R

script.arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
# Repository paths are defined by R/project_paths.R.

plot.style <- Sys.getenv("CRC_PLOT_STYLE", unset = "violin")
if (!plot.style %in% c("violin", "logratio"))
  stop("CRC_PLOT_STYLE must be 'violin' or 'logratio'.")
reuse.summary <- plot.style == "logratio" ||
  identical(Sys.getenv("CRC_REUSE_SUMMARY", unset = "false"), "true")
required <- c("ggplot2", "forestplot", "gridExtra", "ggplotify")
if (!reuse.summary) required <- c("MetaQCAT", required)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

processed <- new.env(parent = emptyenv())
load(file.path(project.dir, "res", "CRC", "wild_score_bootstrap",
               "CRC_processed.RData"),
     envir = processed)
latest <- new.env(parent = emptyenv())
load(file.path(project.dir, "res", "CRC", "wild_score_bootstrap", "CRC_MetaQCAT_bootstrap_res.RData"),
     envir = latest)
if (!reuse.summary) {
  meta.sub <- new.env(parent = globalenv())
  sys.source(file.path(project.dir, "code", "MetaSub.R"), envir = meta.sub)
  if (!exists("MetaSub_latest", envir = meta.sub, inherits = FALSE))
    stop("MetaSub_latest() is missing from code/MetaSub.R")
}

families <- c("543 Enterobacteriaceae", "186803 Lachnospiraceae",
              "186804 Peptostreptococcaceae", "31979 Clostridiaceae")
requested.families <- commandArgs(trailingOnly = TRUE)
if (length(requested.families)) families <- requested.families
reference.override <- Sys.getenv("CRC_PLOT_REFERENCE_TAXON", unset = "")
if (!nzchar(reference.override)) reference.override <- NULL
if (!is.null(reference.override) && length(families) != 1L)
  stop("Select exactly one family when overriding the reference taxon.")
study.names <- names(processed$counts)
if (!identical(study.names, names(processed$disease)) ||
    !identical(study.names, c("AUS", "CHI", "FRA", "GER", "USA"))) {
  stop("Processed CRC count/disease study lists are not aligned")
}
for (study in study.names) {
  if (!identical(rownames(processed$counts[[study]]),
                 rownames(processed$disease[[study]]))) {
    stop("Processed CRC sample rows are not aligned in ", study)
  }
}
if (!all(families %in% colnames(latest$FE_CRC$lineage.pval)) ||
    !all(families %in% colnames(latest$RE_CRC$lineage.pval))) {
  stop("A selected family is absent from the latest CRC result")
}

result.dir <- file.path(project.dir, "outcome", "CRC", "selected_O_lineages")
plot.dir <- file.path(result.dir, "plots")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

group.info <- do.call(rbind, lapply(study.names, function(study) {
  status <- processed$disease[[study]][, 1L]
  if (anyNA(status) || !all(status %in% c(0, 1)))
    stop("Disease indicator must be binary in ", study)
  data.frame(study = study, Group = ifelse(status == 1, "CRC", "CTR"),
             stringsAsFactors = FALSE)
}))
group.info$study <- factor(group.info$study, levels = study.names)
group.info$Group <- factor(group.info$Group, levels = c("CRC", "CTR"))

plot.label <- function(taxon) {
  if (grepl("^NA Enterobacteriaceae gen\\. \\[C ", taxon)) {
    return("Escherichia/Shigella")
  }
  if (grepl("^NA .* gen\\. incertae sedis$", taxon)) return("Incertae sedis")
  if (grepl("^NA unknown ", taxon)) return("Unclassified")
  sub("^[0-9]+ ", "", taxon)
}

make_forest <- function(effects, taxon) {
  rows <- data.frame(study.id = study.names, mean = NA_real_,
                     lower = NA_real_, upper = NA_real_)
  observed <- effects[effects$taxon == taxon & effects$study != "Overall (FE)", ]
  if (nrow(observed)) {
    at <- match(observed$study, study.names)
    if (anyNA(at) || anyDuplicated(at)) stop("Unaligned study effects for ", taxon)
    rows$mean[at] <- observed$estimate
    rows$lower[at] <- observed$ci_lower
    rows$upper[at] <- observed$ci_upper
  }
  pooled <- effects[effects$taxon == taxon & effects$study == "Overall (FE)", ]
  if (nrow(pooled) != 1L) stop("Missing pooled FE effect for ", taxon)
  # forestplot uses the first row to establish its coordinate system; the
  # absent-study sentinel matches the earlier CRC summary-plot convention.
  if (is.na(rows$mean[[1L]])) rows$mean[[1L]] <- -100
  forest <- forestplot::forestplot(
    rows, mean = mean, lower = lower, upper = upper, labeltext = study.id,
    txt_gp = forestplot::fpTxtGp(xlab = grid::gpar(cex = 0.75)),
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
  if (length(at) != 1L) stop("Expected one visible CRC/CTR legend")
  built$grobs[[at]]
}

plot.family <- function(family, summary, file.stem) {
  child.counts <- summary$counts
  rarefied <- summary$rarefied_counts
  if (!identical(names(child.counts), study.names) ||
      !identical(names(rarefied), study.names)) {
    stop("MetaSub study order differs from processed CRC data for ", family)
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
  expected.ref <- if (is.null(reference.override)) {
    latest$FE_CRC$lineage.coef[[family]]$reference.taxon
  } else {
    reference.override
  }
  if (!identical(reference, expected.ref) || !reference %in% taxa)
    stop("MetaSub reference taxon differs from the requested reference for ", family)
  taxon.order <- taxa[order(colSums(Y), decreasing = TRUE)]
  taxon.order <- c(taxon.order[taxon.order != reference], reference)
  if (anyDuplicated(vapply(taxon.order, plot.label, character(1))))
    stop("Plot genus labels are not unique for ", family)

  positive.effects <- summary$positive$estimates
  zero.effects <- summary$zero$estimates
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
  # Spread a large family evenly across pages of at most five taxa; this
  # avoids a nearly empty final page when the number is not a multiple of five.
  # Keep one legend at the last nonempty chart of each component on each page.
  n.taxa <- length(taxon.order)
  n.pages <- ceiling(n.taxa / 5L)
  page.sizes <- rep(n.taxa %/% n.pages, n.pages)
  if (n.taxa %% n.pages)
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

  palette <- c(CRC = "#CD3333", CTR = "steelblue")
  quartiles <- function(x) c(ymin = stats::quantile(x, .25, na.rm = TRUE),
                             ymax = stats::quantile(x, .75, na.rm = TRUE))
  # Keep each chart in the same-sized cell. A legend gets its own narrow
  # sidecar beside the last nonempty chart, never reducing the chart width.
  chart.width <- 4.0
  legend.width <- 0.9
  chart.inset <- 0.1
  column.gap <- 0.22
  axis.title.width <- 0.62
  heading.height <- 0.55
  chart.height <- 3.75
  forest.height <- 2.75
  row.heights <- grid::unit(c(heading.height, chart.height, forest.height,
                              chart.height, forest.height), "in")
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
        if (any(pooled$numerator <= 0 | pooled$ref.count <= 0))
          stop("A pooled taxon/reference count is zero for ", family,
               " / ", taxon)
        pooled$value <- log(pooled$numerator / pooled$ref.count)
        pooled$x <- as.numeric(pooled$study) +
          ifelse(pooled$Group == "CRC", -.18, .18)
        crc <- pooled[pooled$Group == "CRC", c("study", "value")]
        ctr <- pooled[pooled$Group == "CTR", c("study", "value")]
        paired <- merge(crc, ctr, by = "study", suffixes = c(".CRC", ".CTR"))
        links <- data.frame(x = as.numeric(paired$study) - .18,
                            xend = as.numeric(paired$study) + .18,
                            y = paired$value.CRC, yend = paired$value.CTR)
        dots <- d[is.finite(d$value), , drop = FALSE]
        dots$x <- as.numeric(dots$study) +
          ifelse(dots$Group == "CRC", -.18, .18)
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
            linewidth = .65, colour = "#444444") +
          ggplot2::geom_point(
            data = pooled,
            ggplot2::aes(x = x, y = value, fill = Group),
            shape = 21, size = 3.2, stroke = .7, colour = "black") +
          ggplot2::scale_x_continuous(
            breaks = seq_along(study.names), labels = study.names,
            limits = c(.5, length(study.names) + .5)) +
          ggplot2::scale_size_continuous(range = c(.45, 2.2), guide = "none") +
          ggplot2::scale_colour_manual(values = palette, guide = "none") +
          ggplot2::scale_fill_manual(values = palette) +
          ggplot2::theme_bw() +
          ggplot2::labs(y = NULL, x = NULL) +
          ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.key.size = grid::unit(.75, "lines"),
            aspect.ratio = 1,
            plot.margin = ggplot2::margin(.1, 0, .1, .6, "cm")
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
          ggplot2::stat_summary(fun = quartiles, geom = "line", linewidth = .5,
                                position = ggplot2::position_dodge(width = .9),
                                color = "white", alpha = .5,
                                show.legend = FALSE, na.rm = TRUE) +
          ggplot2::stat_summary(fun = "mean", geom = "point", shape = 2,
                                position = ggplot2::position_dodge(width = .9),
                                show.legend = FALSE, na.rm = TRUE) +
          ggplot2::scale_x_discrete(drop = FALSE) +
          ggplot2::scale_fill_manual(values = palette) +
          ggplot2::theme_bw() +
          ggplot2::labs(y = NULL, x = NULL) +
          ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.key.size = grid::unit(.75, "lines"),
            aspect.ratio = 1,
            plot.margin = ggplot2::margin(.1, 0, .1, .6, "cm")
          )
      }
      if (taxon %in% positive.legend)
        pos.legend.grob <- extract.legend(violin)
      violin <- violin + ggplot2::theme(legend.position = "none")
      violin <- gridExtra::arrangeGrob(
        grid::nullGrob(), violin, nrow = 1,
        widths = grid::unit(c(chart.inset, chart.width - chart.inset), "in"))
      pos.forest <- make_forest(positive.effects, taxon)
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
        ggplot2::theme_bw() +
        ggplot2::labs(y = NULL, x = NULL) +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.key.size = grid::unit(.75, "lines"),
          aspect.ratio = 1,
          plot.margin = ggplot2::margin(.1, 0, .1, .2, "cm")
        )
      if (taxon %in% zero.legend)
        zero.legend.grob <- extract.legend(bars)
      bars <- bars + ggplot2::theme(legend.position = "none")
      bars <- gridExtra::arrangeGrob(
        grid::nullGrob(), bars, nrow = 1,
        widths = grid::unit(c(chart.inset, chart.width - chart.inset), "in"))
      zero.forest <- make_forest(zero.effects, taxon)
    } else {
      bars <- empty.panel()
      zero.forest <- empty.panel()
    }
    title <- paste0(plot.label(taxon),
                    if (identical(taxon, reference)) " (Ref.)" else "")
    chart.column <- gridExtra::arrangeGrob(
      grid::textGrob(title, x = grid::unit(.55, "npc"),
                     gp = grid::gpar(fontsize = 18, fontface = "bold")),
      violin, pos.forest, bars, zero.forest,
      ncol = 1, heights = row.heights
    )
    if (taxon %in% positive.legend || taxon %in% zero.legend) {
      legend.column <- gridExtra::arrangeGrob(
        grid::nullGrob(), pos.legend.grob, grid::nullGrob(),
        zero.legend.grob, grid::nullGrob(),
        ncol = 1, heights = row.heights
      )
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
  pages <- lapply(page.indices, function(index) {
    axis.titles <- gridExtra::arrangeGrob(
      grid::nullGrob(),
      grid::textGrob(positive.axis.title, rot = 90,
                     gp = grid::gpar(fontsize = 16, fontface = "bold")),
      grid::textGrob("Proportion of zero counts", rot = 90,
                     gp = grid::gpar(fontsize = 16, fontface = "bold")),
      ncol = 1,
      heights = grid::unit(c(heading.height, chart.height + forest.height,
                             chart.height + forest.height), "in"))
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
  figure.height <- sum(c(heading.height, chart.height, forest.height,
                         chart.height, forest.height))
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
  observed.taxa <- unique(unlist(lapply(processed$counts, colnames),
                                 use.names = FALSE))
  genus <- intersect(unique(processed$taxonomy[
    processed$taxonomy[, 5L] == family, 6L]), observed.taxa)
  if (length(genus) < 3L)
    stop("Selected family has fewer than three genera: ", family)
  file.stem <- sub("^[0-9]+ ", "", family)
  if (!is.null(reference.override)) {
    reference.label <- gsub("[^[:alnum:]]+", "_",
                            sub("^[0-9]+ ", "", reference.override))
    file.stem <- paste0(file.stem, "_ref_", reference.label)
  }
  summary.path <- file.path(result.dir, paste0(file.stem, "_summary.rds"))
  if (reuse.summary) {
    if (!file.exists(summary.path))
      stop("Run the original summary before reusing it: ", summary.path)
    summary <- readRDS(summary.path)
    if (!identical(summary$lineage, family) ||
        !identical(summary$min.prev, 0.1) ||
        !identical(summary$min.depth, 0))
      stop("Saved summary does not match the current CRC setting; ",
           "regenerate with min.prev = 0.1: ", family)
  } else {
    summary <- meta.sub$MetaSub_latest(
      counts = processed$counts, X = processed$disease, X.index = 1L,
      taxonomy = processed$taxonomy, lineage = family,
      cluster.id = processed$subject_ids, min.depth = 0,
      min.prev = 0.1, seed = 1234L,
      reference.taxon = reference.override
    )
    saveRDS(summary, summary.path, compress = "xz")
  }
  output.stem <- paste0(file.stem,
                        if (plot.style == "logratio") "_logratio" else "",
                        "_summary")
  plot.family(family, summary, output.stem)
}
