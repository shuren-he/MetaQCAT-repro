#!/usr/bin/env Rscript
# Rscript run_simulation.R /path/config.R smoke /absolute/output
# Run principal/pipeline/power only deliberately: each requests 5000 datasets
# per setting and 1999 bootstrap draws. Nothing runs during R CMD check.
suppressPackageStartupMessages(library(MetaQCAT))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
  stop("Usage: Rscript run_simulation.R config.R smoke|principal|pipeline|power output_dir")
source(args[1L])
cfg <- validation_config(args[2L])
output.dir <- args[3L]
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
settings <- expand.grid(structure = c("independent", "repeated"),
  covariates = c("independent", "associated"),
  size = names(cfg$size.multipliers), scenario = cfg$scenarios,
  signal = cfg$signals, stringsAsFactors = FALSE)
settings$key <- apply(settings, 1L, paste, collapse = "_")
set.seed(cfg$seed, kind = "L'Ecuyer-CMRG")
seeds <- matrix(sample.int(.Machine$integer.max,
  2L * nrow(settings) * cfg$n.rep, replace = FALSE), ncol = 2L,
  dimnames = list(NULL, c("data", "bootstrap")))

# Saturated stratum models represent both count means and zero probabilities.
# Conditional on stratum, disease is independent of latent composition under
# the null. Cluster size is independent of disease/latent composition given
# stratum. Repeated visits share the subject's Dirichlet composition.
generate_input <- function(setting, seed) {
  set.seed(seed)
  counts <- X <- ids <- vector("list", length(cfg$subject.counts))
  names(counts) <- names(X) <- names(ids) <- paste0("Study", seq_along(counts))
  n.subject <- cfg$subject.counts * cfg$size.multipliers[[setting$size]]
  target <- sample.int(cfg$taxa, 1L)
  direction <- sample(c(-1, 1), 1L)
  flip.study <- sample.int(length(counts), 1L)
  for (s in seq_along(counts)) {
    g <- n.subject[s]
    stratum <- sample(rep(0:1, length.out = g))
    propensity <- if (setting$covariates == "associated")
      cfg$disease.probability.associated else cfg$disease.probability.independent
    disease <- rbinom(g, 1L, propensity[stratum + 1L])
    visits <- if (setting$structure == "repeated")
      sample(cfg$repeat.sizes, g, replace = TRUE) else rep(1L, g)
    subject <- rep(seq_len(g), visits)
    mean0 <- c(.4, .3, .2, .1)
    shift <- c(.25, -.1, .1, -.25) * (1 + .1 * s)
    latent <- matrix(0, g, cfg$taxa)
    effect.direction <- direction * ifelse(setting$scenario == "one_flipped" &&
                                              s == flip.study, -1, 1)
    for (i in seq_len(g)) {
      p <- mean0 * exp(stratum[i] * shift)
      p <- p / sum(p)
      p <- rgamma(cfg$taxa, shape = cfg$dirichlet.concentration * p)
      p <- p / sum(p)
      if (setting$signal == "count")
        p[target] <- p[target] * exp(cfg$signal.log.multiplier *
                                     (disease[i] == as.integer(effect.direction > 0)))
      latent[i, ] <- p / sum(p)
    }
    depth <- if (cfg$pipeline == "fixed") rep(cfg$fixed.depth, length(subject)) else
      pmax(1L, rpois(length(subject), cfg$fixed.depth * (1 + .2 * s)))
    y <- t(vapply(seq_along(subject), function(v) {
      out <- as.vector(rmultinom(1L, depth[v], latent[subject[v], ]))
      if (setting$signal == "zero" &&
          disease[subject[v]] == as.integer(effect.direction > 0)) {
        # A zero-driven perturbation, not a pure zero-only alternative.
        out[target] <- out[target] + rbinom(1L, 5L, cfg$signal.add.probability)
        out <- as.vector(rmultinom(1L, depth[v], out))
      }
      out
    }, numeric(cfg$taxa)))
    dimnames(y) <- list(paste0("sample", seq_along(subject)),
                         paste0("taxon", seq_len(cfg$taxa)))
    x <- cbind(disease = disease[subject], stratum = stratum[subject])
    rownames(x) <- rownames(y)
    counts[[s]] <- y; X[[s]] <- x
    ids[[s]] <- factor(paste0("subject", subject))
  }
  list(counts = counts, X = X, cluster.id = ids,
       n.subject = sum(n.subject), n.sample = sum(vapply(counts, nrow, integer(1))))
}

# Fixed benchmark: fit each component once, then reuse its cached summary
# for FE and RE with a common keyed multiplier plan. No rarefaction or extra
# filtering. Unsupported eligibility is recorded as failure rather than being
# silently interpreted as the fixed-feature benchmark.
fit_fixed <- function(input, seed) {
  ids <- MetaQCAT:::.complete_cluster_ids(input$counts, input$cluster.id)
  weights <- MetaQCAT:::.make_cluster_weights(ids, cfg$n.boot, seed)
  x <- lapply(input$X, function(z) cbind(intercept = 1, z))
  out <- setNames(rep(NA_real_, length(cfg$methods)), cfg$methods)
  reasons <- setNames(rep(NA_character_, length(out)), names(out))
  for (component in c("count", "zero")) {
    summary <- tryCatch({
      if (component == "count") {
        if (any(vapply(input$counts, function(y) any(colSums(y) == 0), logical(1))))
          stop("Fixed benchmark has an unmeasured all-zero feature")
        MetaQCAT:::.Score_test_QCAT_stat(input$counts, x, 2L,
          prev.filter = 0, cluster.list = ids, meta.method = "RE-MetaQCAT",
          reference.taxon = "taxon4")
      } else {
        z <- lapply(input$counts, function(y) (y == 0) * 1)
        if (any(vapply(z, function(y) any(colMeans(y) %in% c(0, 1)), logical(1))))
          stop("Fixed benchmark has an invariant zero feature")
        MetaQCAT:::.Score_test_zero_stat(z, x, 2L, prev.filter = 0,
          cluster.list = ids, meta.method = "RE-MetaQCAT*")
      }
    }, error = identity)
    methods <- paste0(c("FE-MetaQCAT", "RE-MetaQCAT"),
                       if (component == "zero") "*" else "")
    for (m in methods) {
      result <- tryCatch({
        if (inherits(summary, "error")) stop(conditionMessage(summary))
        dimension <- if (component == "count") summary$n.par.interest.beta else
          summary$n.par.interest.alpha
        b <- MetaQCAT:::score_bootstrap_meta(summary,
          MetaQCAT:::.weights_for_summary(summary, weights),
          summary$col.index.list, dimension, 1L, m, component)
        if (length(b$bootstrap.statistic) != cfg$n.boot ||
            any(!is.finite(b$bootstrap.statistic))) stop("Incomplete/nonfinite draws")
        (1 + sum(b$bootstrap.statistic >= b$observed.statistic)) / (cfg$n.boot + 1)
      }, error = identity)
      if (inherits(result, "error")) reasons[m] <- conditionMessage(result) else out[m] <- result
    }
  }
  for (model in c("FE", "RE")) {
    components <- paste0(model, c("-MetaQCAT", "-MetaQCAT*"))
    m <- paste0(model, "-MetaQCAT-O")
    if (all(is.finite(out[components]))) out[m] <- MetaQCAT:::ACAT(out[components]) else
      reasons[m] <- paste(na.omit(reasons[components]), collapse = "; ")
  }
  list(p = out, reason = reasons)
}

records <- vector("list", nrow(settings) * cfg$n.rep)
for (s in seq_len(nrow(settings))) for (r in seq_len(cfg$n.rep)) {
  index <- (s - 1L) * cfg$n.rep + r
  input <- generate_input(settings[s, ], seeds[index, "data"])
  start <- proc.time()[["elapsed"]]
  if (cfg$pipeline == "fixed") {
    fit <- fit_fixed(input, seeds[index, "bootstrap"])
  } else {
    p <- reason <- setNames(rep(NA_real_, length(cfg$methods)), cfg$methods)
    reason <- setNames(rep(NA_character_, length(p)), names(p))
    # Exactly one generated dataset, with the same seed/rarefaction for every
    # method. Full preprocessing is repeated only for a new dataset.
    for (m in cfg$methods) {
      value <- tryCatch(suppressWarnings(MetaQCAT(input$counts, input$X, 1L,
        cluster.id = input$cluster.id, method = m, min.depth = cfg$min.depth,
        min.prev = cfg$min.prev, n.boot = cfg$n.boot, seed = seeds[index, "bootstrap"],
        n.threads = cfg$n.threads)), error = identity)
      field <- if (grepl("-O$", m)) "pvalue.comb" else
        if (grepl("\\*$", m)) "pvalue.pa" else "pvalue.comp"
      if (inherits(value, "error")) reason[m] <- conditionMessage(value) else {
        p[m] <- value[[field]][["Boot"]]
        if (!is.finite(p[m])) reason[m] <- value$error
      }
    }
    fit <- list(p = p, reason = reason)
  }
  records[[index]] <- data.frame(setting = settings$key[s], replicate = r,
    method = cfg$methods, p.value = unname(fit$p[cfg$methods]),
    failure.reason = unname(fit$reason[cfg$methods]),
    runtime.dataset.seconds = proc.time()[["elapsed"]] - start,
    n.subject = input$n.subject, n.sample = input$n.sample,
    n.boot.requested = cfg$n.boot,
    n.boot.completed = ifelse(is.finite(fit$p[cfg$methods]), cfg$n.boot, 0L),
    data.seed = seeds[index, "data"], bootstrap.seed = seeds[index, "bootstrap"])
  message(settings$key[s], ": ", r, "/", cfg$n.rep)
}
raw <- do.call(rbind, records)
raw$failed <- !is.finite(raw$p.value)
raw$nonestimable <- raw$failed & grepl(
  "rank|singular|condition|information|invariant zero feature|all-zero feature|No studies",
  ifelse(is.na(raw$failure.reason), "", raw$failure.reason), ignore.case = TRUE)
summary <- do.call(rbind, lapply(split(raw, list(raw$setting, raw$method), drop = TRUE), function(x) {
  do.call(rbind, lapply(cfg$alpha, function(alpha) {
    valid <- !x$failed; n <- sum(valid); k <- sum(x$p.value[valid] < alpha)
    ci <- if (n) binom.test(k, n)$conf.int else c(NA_real_, NA_real_)
    data.frame(setting = x$setting[1L], method = x$method[1L], alpha = alpha,
      attempted = nrow(x), successful = n, failures = sum(x$failed),
      nonestimable = sum(x$nonestimable),
      other.failures = sum(x$failed & !x$nonestimable),
      failure.rate = mean(x$failed), rejections = k,
      rejection.rate = if (n) k / n else NA_real_,
      mc.se = if (n) sqrt((k/n)*(1-k/n)/n) else NA_real_,
      mc.lower = ci[1L], mc.upper = ci[2L],
      rejection.rate.all.attempted = k / nrow(x),
      mean.dataset.runtime.seconds = mean(x$runtime.dataset.seconds),
      n.boot = cfg$n.boot, n.perm = 0L,
      status = if (cfg$n.rep < 100L) "execution_smoke_only" else "requires_statistical_review")
  }))
}))
rownames(summary) <- NULL
result <- list(configuration = cfg, settings = settings, seeds = seeds,
  package.version = as.character(packageVersion("MetaQCAT")),
  code.revision = unname(tools::md5sum(system.file("validation", "source_manifest.csv",
                                                 package = "MetaQCAT"))),
  installed.package.path = find.package("MetaQCAT"), session.info = sessionInfo(),
  source.md5 = {
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    tools::md5sum(c(args[1L], script))
  },
  raw = raw, summary = summary,
  note = paste("Binomial intervals use successful fits; all attempted fits and",
    "failures are reported separately. Runtime is for a paired six-method dataset,",
    "not a separate per-method fit cost. No FDR-control claim; no automatic",
    "pass/fail assessment based only on whether an interval contains alpha."))
saveRDS(result, file.path(output.dir, "validation.rds"))
write.csv(raw, file.path(output.dir, "replicates.csv"), row.names = FALSE)
write.csv(summary, file.path(output.dir, "summary.csv"), row.names = FALSE)
write.csv(subset(raw, failed), file.path(output.dir, "failures.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
