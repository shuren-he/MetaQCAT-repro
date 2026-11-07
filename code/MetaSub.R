# Taxon-level one-step score summaries for the paper figures.
# Use the same count/zero score backend as the current MetaQCAT package.
# No separate permutation or C++ implementation is loaded here.
if (!exists("project.dir", inherits = TRUE))
  source(file.path("R", "project_paths.R"))

MetaSub_latest <- function(counts, X, X.index, taxonomy, lineage,
                           cluster.id = NULL, min.depth = 0,
                           min.prev = 0.1, seed = 1234L,
                           output.dir = NULL, reference.taxon = NULL) {
  if (!requireNamespace("MetaQCAT", quietly = TRUE))
    stop("Install the current MetaQCAT package before using MetaSub_latest().")
  if (utils::packageVersion("MetaQCAT") < "0.2.1")
    stop("MetaSub_latest() requires MetaQCAT 0.2.1 or later.")
  if (length(X.index) != 1L)
    stop("MetaSub_latest() currently plots one tested covariate at a time.")
  if (!is.character(lineage) || length(lineage) != 1L ||
      is.na(lineage) || !nzchar(lineage))
    stop("lineage must be one nonempty family name.")
  if (!is.numeric(min.depth) || length(min.depth) != 1L ||
      !is.finite(min.depth) || min.depth < 0)
    stop("min.depth must be one nonnegative number.")

  tax <- as.matrix(taxonomy)
  if (ncol(tax) == 2L) {
    colnames(tax) <- c("Rank1", "Rank2")
  } else if (ncol(tax) >= 6L) {
    tax <- tax[, 5:6, drop = FALSE]
    colnames(tax) <- c("Rank1", "Rank2")
  } else {
    stop("taxonomy must contain family and genus ranks in two columns or ranks 5:6.")
  }
  prepare <- getFromNamespace(".prepare_inputs", "MetaQCAT")
  rarefy <- getFromNamespace(".rarefy_counts", "MetaQCAT")
  with.seed <- getFromNamespace(".with_seed", "MetaQCAT")
  d <- prepare(counts, X, X.index, cluster.id = cluster.id, taxonomy = tax)
  for (i in seq_along(d$counts)) {
    keep <- rowSums(d$counts[[i]]) > min.depth
    d$counts[[i]] <- d$counts[[i]][keep, , drop = FALSE]
    d$X[[i]] <- d$X[[i]][keep, , drop = FALSE]
    if (!is.null(d$cluster.id) && !is.null(d$cluster.id[[i]]))
      d$cluster.id[[i]] <- factor(d$cluster.id[[i]][keep])
  }
  keep.study <- vapply(d$counts, nrow, integer(1)) > 0L
  if (!any(keep.study)) stop("No studies remain after depth filtering.")
  d$counts <- d$counts[keep.study]
  d$X <- d$X[keep.study]
  if (!is.null(d$cluster.id)) d$cluster.id <- d$cluster.id[keep.study]

  family.names <- unique(d$taxonomy[, "Rank1"])
  family.key <- lineage
  if (!family.key %in% family.names) {
    plain <- sub("^f__", "", lineage)
    matched <- family.names[sub("^[0-9]+[[:space:]]+", "", family.names) == plain]
    if (length(matched) != 1L)
      stop("Family not found uniquely in taxonomy: ", lineage)
    family.key <- matched
  }
  rows <- which(d$taxonomy[, "Rank1"] == family.key)
  child <- d$taxonomy[rows, "Rank2"]
  if (length(unique(child)) < 2L)
    stop("The selected family has fewer than two genus-level taxa.")
  aggregate.family <- function(m) {
    out <- t(rowsum(t(m[, rows, drop = FALSE]), child, reorder = FALSE))
    rownames(out) <- rownames(m)
    storage.mode(out) <- "double"
    out
  }
  Y <- lapply(d$counts, aggregate.family)
  # Match the public call: one RNG stream across studies, before aggregation.
  Y.rare <- lapply(with.seed(
    as.integer((as.double(seed) + 1) %% .Machine$integer.max),
    lapply(d$counts, rarefy)), aggregate.family)
  design <- lapply(d$X, function(x) cbind(`(Intercept)` = 1, x))
  tested <- d$X.index + 1L
  positive.stat <- if (is.null(reference.taxon)) {
    MetaQCAT::Score_test_QCAT_stat(
      counts = Y, X = design, X.index = tested, min.prev = min.prev,
      cluster.id = d$cluster.id, method = "FE-MetaQCAT")
  } else {
    # The public score wrapper does not expose its reference-taxon argument.
    getFromNamespace(".Score_test_QCAT_stat", "MetaQCAT")(
      Y.list = Y, X.list = design, X.par.index = tested,
      prev.filter = min.prev, cluster.list = d$cluster.id,
      meta.method = "FE-MetaQCAT", reference.taxon = reference.taxon)
  }
  zero.indicator <- lapply(Y.rare, function(y) (y == 0) * 1)
  zero.stat <- MetaQCAT::Score_test_zero_stat(
    zeros = zero.indicator, X = design, X.index = tested,
    min.prev = min.prev, cluster.id = d$cluster.id,
    method = "FE-MetaQCAT*")

  reference.taxon <- positive.stat$reference.taxon
  taxon.order <- c(setdiff(colnames(Y[[1L]]), reference.taxon), reference.taxon)
  Y <- lapply(Y, function(y) y[, taxon.order, drop = FALSE])
  Y.rare <- lapply(Y.rare, function(y) y[, taxon.order, drop = FALSE])

  summarize.part <- function(stat, part) {
    is.positive <- identical(part, "positive")
    local.y <- if (is.positive) stat$Y.list else stat$Y.bin.list
    local.est <- if (is.positive) stat$est.beta.lst else stat$est.alpha.lst
    pooled.est <- if (is.positive) stat$beta.hat else stat$alpha.hat
    pooled.var <- if (is.positive) stat$var.beta.hat else stat$var.alpha.hat
    if (length(local.y) != length(local.est) ||
        length(local.y) != length(stat$cov.lst))
      stop("Backend study estimates are not aligned for ", part, ".")
    study.names <- names(local.y)
    if (is.null(study.names) || anyNA(study.names) ||
        !all(study.names %in% names(Y)))
      stop("Backend did not retain study names for ", part, ".")
    pooled.taxa <- colnames(pooled.est)
    if (is.null(pooled.taxa) || anyNA(pooled.taxa))
      stop("Backend pooled taxa are unnamed for ", part, ".")
    effect.rows <- vector("list", length(local.y) + 1L)
    for (i in seq_along(local.y)) {
      taxa <- colnames(local.y[[i]])
      if (is.positive) taxa <- taxa[-length(taxa)]
      estimate <- as.numeric(local.est[[i]])
      variance <- diag(as.matrix(stat$cov.lst[[i]]))
      if (length(taxa) != length(estimate) ||
          length(taxa) != length(variance) ||
          any(!is.finite(estimate)) || any(!is.finite(variance)) ||
          any(variance < -1e-10))
        stop("Invalid study estimate/covariance for ", part,
             " in ", study.names[[i]], ".")
      effect.rows[[i]] <- data.frame(
        study = study.names[[i]], taxon = taxa,
        estimate = estimate, se = sqrt(pmax(variance, 0)),
        stringsAsFactors = FALSE)
    }
    pooled.estimate <- as.numeric(pooled.est)
    pooled.variance <- as.numeric(pooled.var)
    if (length(pooled.taxa) != length(pooled.estimate) ||
        any(!is.finite(pooled.estimate)) ||
        any(!is.finite(pooled.variance)) ||
        any(pooled.variance < -1e-10))
      stop("Invalid pooled estimate/covariance for ", part, ".")
    effect.rows[[length(effect.rows)]] <- data.frame(
      study = "Overall (FE)", taxon = pooled.taxa,
      estimate = pooled.estimate, se = sqrt(pmax(pooled.variance, 0)),
      stringsAsFactors = FALSE)
    effects <- do.call(rbind, effect.rows)
    rownames(effects) <- NULL
    effects$ci_lower <- effects$estimate - 1.96 * effects$se
    effects$ci_upper <- effects$estimate + 1.96 * effects$se

    # Reconstruct the full pooled covariance from the study precisions. Its
    # diagonal is fixed to the package-returned marginal variances so the
    # forest intervals exactly follow the current package result.
    information <- matrix(0, length(pooled.taxa), length(pooled.taxa))
    for (i in seq_along(local.y)) {
      taxa <- colnames(local.y[[i]])
      if (is.positive) taxa <- taxa[-length(taxa)]
      index <- match(taxa, pooled.taxa)
      if (anyNA(index)) stop("Local taxa absent from pooled order.")
      information[index, index] <- information[index, index] +
        MASS::ginv(as.matrix(stat$cov.lst[[i]]))
    }
    cov.mu <- MASS::ginv(information)
    diag(cov.mu) <- pooled.variance
    dimnames(cov.mu) <- list(pooled.taxa, pooled.taxa)
    re.study <- lapply(seq_along(local.y), function(i) {
      tryCatch({
        if (is.positive) {
          re <- getFromNamespace("re_qcat_components", "MetaQCAT")(
            as.numeric(local.est[[i]]), as.matrix(stat$cov.lst[[i]]), 1L)
          re$statistic <- max(0, re$U_theta)^2 / re$V_theta
          re
        } else {
          V <- as.matrix(stat$cov.lst[[i]])
          precision <- MASS::ginv(V)
          z <- precision %*% as.numeric(local.est[[i]])
          U <- 0.5 * (sum(z * z) - sum(diag(precision)))
          V.theta <- 0.5 * sum(precision * precision)
          if (!is.finite(V.theta) || V.theta <= 0)
            stop("Degenerate zero-part RE information.")
          list(U_theta = U, V_theta = V.theta,
               statistic = max(0, U)^2 / V.theta)
        }
      }, error = function(e) list(error = conditionMessage(e)))
    })
    names(re.study) <- study.names
    if (any(vapply(re.study, function(z) !is.null(z$error), logical(1)))) {
      random.effect <- list(error = paste(vapply(re.study, function(z)
        if (is.null(z$error)) "" else z$error, character(1)), collapse = "; "),
        by.study = re.study)
    } else {
      U <- sum(vapply(re.study, `[[`, numeric(1), "U_theta"))
      V.theta <- sum(vapply(re.study, `[[`, numeric(1), "V_theta"))
      random.effect <- list(U_theta = U, V_theta = V.theta,
        statistic = max(0, U)^2 / V.theta,
        by.study = re.study)
    }
    list(study_names = study.names, taxa = taxon.order,
         estimates = effects, score.stat.meta = stat$score.stat.meta,
         reference.taxon = if (is.positive) reference.taxon else NULL,
         random_effect = random.effect, pooled_covariance = cov.mu,
         backend = stat)
  }

  positive <- summarize.part(positive.stat, "positive")
  zero <- summarize.part(zero.stat, "zero")
  result <- list(counts = Y, rarefied_counts = Y.rare,
                 positive = positive, zero = zero,
                 lineage = lineage, min.prev = min.prev,
                 min.depth = min.depth, seed = seed)
  if (!is.null(output.dir)) {
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    save(Y, file = file.path(output.dir, paste0(lineage, "_count.Rdata")))
    save(Y.rare, file = file.path(output.dir,
                                 paste0(lineage, "_zero_count.Rdata")))
    save.study.id <- match(positive$study_names, names(Y))
    score.beta <- positive.stat$est.beta.lst
    est.cov <- positive.stat$cov.lst
    idx.list <- lapply(positive.stat$Y.list, function(y)
      match(colnames(y)[-ncol(y)], taxon.order))
    save.index.pos <- match(colnames(positive.stat$beta.hat), taxon.order)
    mu.hat <- as.numeric(positive.stat$beta.hat)
    cov.mu <- positive$pooled_covariance
    random.effect <- positive$random_effect
    save(save.study.id, score.beta, est.cov, idx.list, save.index.pos,
         mu.hat, cov.mu, reference.taxon, random.effect,
         file = file.path(output.dir, paste0(lineage, "_pos_meta_info.Rdata")))
    save.study.id <- match(zero$study_names, names(Y))
    score.alpha <- zero.stat$est.alpha.lst
    est.cov.zero <- zero.stat$cov.lst
    idx.list <- lapply(zero.stat$Y.bin.list, function(y)
      match(colnames(y), taxon.order))
    save.index.zero <- match(colnames(zero.stat$alpha.hat), taxon.order)
    mu.hat <- as.numeric(zero.stat$alpha.hat)
    cov.mu <- zero$pooled_covariance
    random.effect <- zero$random_effect
    save(save.study.id, score.alpha, est.cov.zero, idx.list,
         save.index.zero, mu.hat, cov.mu, random.effect,
         file = file.path(output.dir, paste0(lineage, "_zero_meta_info.Rdata")))
  }
  result
}
