# IBD: P/CD/RD cluster-QCAT unions versus FE/RE-MetaQCAT-O bootstrap.
# Run from the repository root after IBD_all_methods.R. No models are fitted.
source(file.path("R", "project_paths.R"))

if (!requireNamespace("UpSetR", quietly = TRUE)) stop("Install UpSetR first.")

result.dir <- Sys.getenv("IBD_QCAT_OUTPUT", unset = file.path(
  project.dir, "res", "IBD", "cluster_qcat"))
q.file <- file.path(result.dir, "IBD_QCAT_cluster_q_val.csv")
q <- read.csv(q.file, row.names = 1L, check.names = FALSE)
baseline.columns <- unlist(lapply(c("RD", "P", "CD"), function(strategy)
  paste0(strategy, "-QCAT", 1:2, " (perm)")))
if (anyDuplicated(rownames(q)) || !all(baseline.columns %in% names(q)))
  stop("Unexpected cluster-QCAT q-value table: ", q.file)
q[baseline.columns] <- lapply(q[baseline.columns], function(x) {
  # read.csv imports a completely missing numeric column as logical.
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  x
})
if (!all(vapply(q[baseline.columns], is.numeric, logical(1))) ||
    any(unlist(q[baseline.columns]) < 0 | unlist(q[baseline.columns]) > 1,
        na.rm = TRUE)) stop("Invalid baseline q-values.")

primary.dir <- file.path(project.dir, "res", "IBD", "wild_score_bootstrap")
fe <- readRDS(file.path(primary.dir, "IBD_primary_FE_MetaQCAT_O.rds"))
re <- readRDS(file.path(primary.dir, "IBD_primary_RE_MetaQCAT_O.rds"))
fe.q <- fe$lineage.qval["Bootstrap-FE-MetaQCAT-O", ]
re.q <- re$lineage.qval["Bootstrap-RE-MetaQCAT-O", ]
lineages <- sort(unique(c(rownames(q), names(fe.q), names(re.q))))

# A strategy detects a lineage if either of its two tests has BH q < 0.05.
# Missing or failed tests are not discoveries; their q-values remain NA.
aligned.q <- matrix(NA_real_, length(lineages), 8L,
  dimnames = list(lineages, c(baseline.columns, "FE-MetaQCAT-O", "RE-MetaQCAT-O")))
aligned.q[match(rownames(q), lineages), baseline.columns] <- as.matrix(q[baseline.columns])
aligned.q[, "FE-MetaQCAT-O"] <- fe.q[match(lineages, names(fe.q))]
aligned.q[, "RE-MetaQCAT-O"] <- re.q[match(lineages, names(re.q))]
sig <- data.frame(Lineage = lineages, check.names = FALSE)
for (strategy in c("RD", "P", "CD")) {
  columns <- paste0(strategy, "-QCAT", 1:2, " (perm)")
  sig[[strategy]] <- as.integer(rowSums(
    !is.na(aligned.q[, columns, drop = FALSE]) &
      aligned.q[, columns, drop = FALSE] < 0.05) > 0L)
}
sig$FE <- as.integer(!is.na(aligned.q[, "FE-MetaQCAT-O"]) &
                       aligned.q[, "FE-MetaQCAT-O"] < 0.05)
sig$RE <- as.integer(!is.na(aligned.q[, "RE-MetaQCAT-O"]) &
                       aligned.q[, "RE-MetaQCAT-O"] < 0.05)
stopifnot(identical(names(sig), c("Lineage", "RD", "P", "CD", "FE", "RE")),
          !anyNA(sig))
write.csv(aligned.q, file.path(result.dir, "IBD_overlap_qvalues.csv"))
write.csv(sig, file.path(result.dir, "IBD_sig_family.csv"), row.names = FALSE)
write.csv(sig, file.path(project.dir, "res", "IBD", "IBD_sig_family.csv"),
          row.names = FALSE)
counts <- data.frame(strategy = names(sig)[-1L],
                     discoveries = as.integer(colSums(sig[-1L])))
write.csv(counts, file.path(result.dir, "IBD_overlap_discovery_counts.csv"),
          row.names = FALSE)
print(counts, row.names = FALSE)

sets <- lapply(sig[-1L], function(x) sig$Lineage[x == 1L])
if (!any(lengths(sets))) stop("No discoveries to display.")
# Identical ordering and plotting settings to CRC_upset_plot.R.
display.order <- names(sets)[order(lengths(sets), seq_along(sets))]
upset.sets <- rev(display.order)
membership <- UpSetR::fromList(sets)
# UpSetR identifies the membership block by its nonconstant edge columns.
# Bracket any all-zero groups with nonempty groups without changing the
# requested display order or adding artificial lineage observations.
positive <- names(sets)[lengths(sets) > 0L]
empty <- names(sets)[lengths(sets) == 0L]
if (length(empty) && length(positive) >= 2L)
  membership <- membership[c(positive[1L], empty, positive[-1L])]
draw <- function() {
  UpSetR::upset(membership,
    sets = upset.sets, keep.order = TRUE, nintersects = 12,
    mb.ratio = c(0.5, 0.5), sets.x.label = "Number of Discoveries",
    order.by = c("freq", "degree"), decreasing = c(TRUE, FALSE),
    point.size = 3, text.scale = c(3, 3, 2, 3, 3, 3))
}
plot.dir <- file.path(project.dir, "outcome", "IBD")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
pdf.file <- file.path(plot.dir, "IBD_overlap_plot.pdf")
png.file <- file.path(plot.dir, "IBD_overlap_plot.png")
grDevices::pdf(pdf.file, width = 11, height = 8)
tryCatch(print(draw(), newpage = FALSE), finally = grDevices::dev.off())
grDevices::png(png.file, width = 11, height = 8, units = "in", res = 300,
               bg = "white")
tryCatch(print(draw(), newpage = FALSE), finally = grDevices::dev.off())
message("Saved cluster-comparator IBD UpSet plot: ", pdf.file)
