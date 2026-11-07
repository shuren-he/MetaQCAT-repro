#!/usr/bin/env Rscript
# Opt-in timing on a real packaged IBD lineage. Run after installing MetaQCAT:
# Rscript benchmark_ibd.R --family=f__Enterobacteriaceae --boot=1000 \
#   --output=/absolute/output/directory
# The script separates null fitting, weight generation, and score evaluation.

suppressPackageStartupMessages(library(MetaQCAT))

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  found <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(found) > 1L) stop("Repeated option: ", name)
  if (!length(found)) return(default)
  sub(paste0("^--", name, "="), "", found)
}
family <- option("family", "f__Enterobacteriaceae")
n.boot <- suppressWarnings(as.integer(option("boot", "1000")))
seed <- suppressWarnings(as.integer(option("seed", "20261101")))
if (is.na(n.boot) || n.boot < 1L || is.na(seed))
  stop("boot must be positive and seed must be an integer.")
output.dir <- option("output", file.path(getwd(), "metaqcat_ibd_benchmark"))
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

data("IBD_data", package = "MetaQCAT")
taxa <- rownames(IBD_data$taxonomy)[
  IBD_data$taxonomy[, "Rank5"] == family
]
if (length(taxa) < 2L) stop("Family has fewer than two taxa: ", family)
align.family <- function(y) {
  aligned <- matrix(0, nrow(y), length(taxa),
                    dimnames = list(rownames(y), taxa))
  present <- intersect(taxa, colnames(y))
  aligned[, present] <- y[, present, drop = FALSE]
  aligned
}
counts <- lapply(IBD_data$counts, align.family)
rarefaction.time <- system.time({
  rarefied <- MetaQCAT:::.with_seed(as.integer((as.double(seed) + 1) %%
    .Machine$integer.max), lapply(IBD_data$counts, MetaQCAT:::.rarefy_counts))
})[["elapsed"]]
zeros <- lapply(rarefied, function(y) (align.family(y) == 0) * 1)
X <- lapply(IBD_data$design, function(x)
  cbind(intercept = 1, x))
ids <- MetaQCAT:::.complete_cluster_ids(counts, IBD_data$subject_ids)
weights <- NULL
weight.time <- system.time({
  weights <- MetaQCAT:::.make_cluster_weights(ids, n.boot, seed)
})[["elapsed"]]

records <- list()
k <- 0L
for (component in c("count", "zero")) {
  stat <- NULL
  fit.time <- system.time({
    stat <- tryCatch({
      if (component == "count") MetaQCAT:::.Score_test_QCAT_stat(
        Y.list = counts, X.list = X, X.par.index = 2L,
        prev.filter = 0.1, cluster.list = ids, meta.method = "RE-MetaQCAT") else
        MetaQCAT:::.Score_test_zero_stat(Y.bin.list = zeros, Z.list = X,
          Z.par.index = 2L, prev.filter = .1, cluster.list = ids,
          meta.method = "RE-MetaQCAT*")
    }, error = function(e) e)
  })[["elapsed"]]
  for (model in c("FE", "RE")) {
  k <- k + 1L
  method <- paste0(model, "-MetaQCAT", if (component == "zero") "*" else "")
  bootstrap.time <- NA_real_
  failure <- NA_character_
  observed <- NA_real_
  if (inherits(stat, "error")) {
    failure <- conditionMessage(stat)
  } else {
    used.weights <- MetaQCAT:::.weights_for_summary(stat, weights)
    result <- NULL
    bootstrap.time <- system.time({
      result <- tryCatch(MetaQCAT:::score_bootstrap_meta(
        stat, used.weights, stat$col.index.list,
        if (component == "count") stat$n.par.interest.beta else
          stat$n.par.interest.alpha, 1L, method, component
      ), error = function(e) e)
    })[["elapsed"]]
    if (inherits(result, "error")) failure <- conditionMessage(result)
    else {
      stopifnot(length(result$bootstrap.statistic) == n.boot)
      observed <- result$observed.statistic
    }
  }
  records[[k]] <- data.frame(
    dataset = "IBD_data", lineage = family, component = component,
    method = method, studies.requested = length(counts),
    samples = sum(vapply(counts, nrow, integer(1))),
    independent.clusters = sum(vapply(ids, function(x)
      length(unique(as.character(x))), integer(1))),
    taxa = length(taxa), null.fit.seconds = fit.time,
    weight.generation.seconds = weight.time,
    rarefaction.seconds = rarefaction.time,
    bootstrap.seconds = bootstrap.time,
    bootstrap.seconds.per.1000 = bootstrap.time * 1000 / n.boot,
    n.boot = n.boot, seed = seed, observed.statistic = observed,
    failed = !is.na(failure), failure.reason = failure,
    stringsAsFactors = FALSE
  )
  }
}
benchmark <- do.call(rbind, records)
result <- list(
  configuration = list(family = family, n.boot = n.boot, seed = seed,
                       X.index = 1L, min.prev = 0.1),
  package.version = as.character(packageVersion("MetaQCAT")),
  code.revision = unname(tools::md5sum(system.file("validation", "source_manifest.csv",
                                                 package = "MetaQCAT"))),
  session.info = sessionInfo(),
  timing = benchmark,
  note = paste(
    "Null fit includes taxon filtering and nuisance-model fitting.",
    "Null-fit timing is shared by FE/RE within a component, not additive.",
    "Weight generation is shared and timed separately.",
    "Bootstrap seconds measure only fixed-fit score evaluation."
  )
)
saveRDS(result, file.path(output.dir, "MetaQCAT_IBD_benchmark.rds"))
write.csv(benchmark, file.path(output.dir, "MetaQCAT_IBD_benchmark.csv"),
          row.names = FALSE)
print(benchmark, row.names = FALSE)
message("Outputs saved to ", normalizePath(output.dir))
