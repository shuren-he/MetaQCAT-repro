# Shared statistical calculations for the six paper simulation entrypoints.
# The generator is the current study-mean-depth experiment. All methods receive
# the same generated counts, labels, retained rows and taxon identities.

r.GDM <- function(X, mod, SeqDepth){
  alpha <- mod$alpha.est
  phi <- mod$phi.est
  if (is.null(alpha) || is.null(phi)) {
    stop("The fitted study model does not contain alpha.est and phi.est.")
  }
  alpha <- as.matrix(alpha)
  N <- nrow(X)
  K <- ncol(alpha)
  if (length(SeqDepth) == 1L) SeqDepth <- rep(SeqDepth, N)
  if (length(SeqDepth) != N) {
    stop("The sequencing-depth vector must have one value per sample.")
  }
  if (ncol(X) != nrow(alpha) || length(phi) != K) {
    stop("The GDM design matrix and fitted parameters have incompatible dimensions.")
  }
  mu.mat <- plogis(X %*% alpha)
  phi.mat <- matrix(phi, nrow = N, ncol = K, byrow = TRUE)
  Z <- matrix(rbeta(N * K, mu.mat * phi.mat,
                    (1 - mu.mat) * phi.mat), N, K)
  P <- matrix(NA_real_, N, K + 1L)
  remaining <- rep(1, N)
  for (j in seq_len(K)) {
    P[, j] <- Z[, j] * remaining
    remaining <- remaining * (1 - Z[, j])
  }
  P[, K + 1L] <- remaining
  P / rowSums(P)
}

scenario_directions <- function(stu.num, type, direct = NULL){
  if (type == "fixed") return(rep(1L, stu.num))
  if (type != "random") stop("type must be 'fixed' or 'random'.")
  if (identical(direct, "both")) {
    ans <- rep(1L, stu.num)
    ans[sample(seq_len(stu.num), 1L)] <- -1L
    return(ans)
  }
  stop("For the heterogeneous power setting, direct must be 'both'.")
}

multinomial_counts <- function(P, depth){
  ans <- t(vapply(seq_len(nrow(P)), function(i) {
    as.integer(rmultinom(1L, size = depth[i], prob = P[i, ])[, 1L])
  }, integer(ncol(P))))
  colnames(ans) <- colnames(P)
  ans
}

renormalize_count_row <- function(x, depth){
  scaled <- x / sum(x) * depth
  ans <- floor(scaled)
  remainder <- as.integer(depth - sum(ans))
  if (remainder > 0L) {
    idx <- order(scaled - ans, decreasing = TRUE)[seq_len(remainder)]
    ans[idx] <- ans[idx] + 1L
  }
  as.integer(ans)
}

simulate_crc_data <- function(sim.cnt, mod_gdm_lst, depth.pool, taxa.names,
                              lst.num, signal.size, signal.method, type,
                              direct = NULL, seed, correlated = FALSE){
  set.seed(seed + sim.cnt)
  stu.num <- length(lst.num)
  n.taxa <- length(taxa.names)
  if (length(mod_gdm_lst) != stu.num || length(depth.pool) != stu.num ||
      any(lengths(depth.pool) == 0L)) {
    stop("GDM models, observed depth pools, and studies must align.")
  }
  if (!all(vapply(mod_gdm_lst, function(mod) {
    length(mod$col.order) == n.taxa
  }, logical(1)))) {
    stop("All fitted GDM study models must contain the same number of taxa.")
  }
  if (!signal.method %in% c("proportion", "count")) {
    stop("signal.method must be 'proportion' or 'count'.")
  }
  if (signal.method == "proportion" &&
      (!is.finite(signal.size) || signal.size < 0)) {
    stop("The proportional increase Delta must be non-negative.")
  }
  if (signal.method == "count" &&
      (!is.finite(signal.size) || signal.size < 0 || signal.size > 1)) {
    stop("The binomial count-signal probability must be between 0 and 1.")
  }
  target.idx <- sample(seq_len(n.taxa), 1L)
  baseline.direction <- sample(c(-1L, 1L), 1L)
  directions <- baseline.direction * scenario_directions(stu.num, type, direct)
  count.trials <- as.integer(round(mean(unlist(depth.pool, use.names = FALSE))))
  dat.lst <- vector("list", stu.num)
  X.lst <- vector("list", stu.num)
  keep.lst <- vector("list", stu.num)
  depth.lst <- vector("list", stu.num)
  cluster.lst <- vector("list", stu.num)
  for (l in seq_len(stu.num)) {
    n <- lst.num[l]
    depth.all = rpois(n, mean(depth.pool[[l]]))
    P.all <- r.GDM(X = cbind(1, rep(0, n)),
      mod = mod_gdm_lst[[l]]$est.para, SeqDepth = depth.all)
    P.all <- P.all[, order(mod_gdm_lst[[l]]$col.order), drop = FALSE]
    colnames(P.all) <- taxa.names
    keep <- which(depth.all >= 1.5 * n.taxa)
    if (length(keep) == 0L) stop(sprintf("Study %d has no samples after the depth filter.", l))
    depth <- depth.all[keep]
    P <- P.all[keep, , drop = FALSE]
    if (correlated) {
      shuffled <- sample.int(length(keep))
      pairs <- split(shuffled, ceiling(seq_along(shuffled) / 2L))
      cluster <- integer(length(keep))
      for (i in seq_along(pairs)) {
        members <- pairs[[i]]
        cluster[members] <- i
        if (length(members) == 2L) {
          shared <- r.GDM(X = matrix(c(1, 0), nrow = 1L),
            mod = mod_gdm_lst[[l]]$est.para, SeqDepth = depth[members[1L]])
          shared <- shared[1L, order(mod_gdm_lst[[l]]$col.order)]
          P[members, ] <- (P[members, , drop = FALSE] +
            matrix(shared, nrow = 2L, ncol = n.taxa, byrow = TRUE)) / 2
        }
      }
      cluster.case <- rep(0L, length(pairs))
      cluster.case[sample.int(length(pairs), ceiling(length(pairs) / 2L))] <- 1L
      X <- cluster.case[cluster]
    } else {
      cluster <- seq_along(keep)
      X <- rep(0L, length(keep))
      X[sample.int(length(keep), ceiling(length(keep) / 2L))] <- 1L
    }
    affected <- if (directions[l] == 1L) X == 1L else X == 0L
    if (directions[l] != 0L && signal.size > 0 && signal.method == "proportion") {
      P[affected, target.idx] <- P[affected, target.idx] * (1 + signal.size)
      P[affected, ] <- P[affected, , drop = FALSE] / rowSums(P[affected, , drop = FALSE])
    }
    Y <- multinomial_counts(P, depth)
    if (directions[l] != 0L && signal.size > 0 && signal.method == "count") {
      affected.idx <- which(affected)
      delta <- rbinom(length(affected.idx), size = count.trials, prob = signal.size)
      for (j in seq_along(affected.idx)) {
        i <- affected.idx[j]
        modified <- Y[i, ]
        modified[target.idx] <- modified[target.idx] + delta[j]
        Y[i, ] <- renormalize_count_row(modified, depth[i])
      }
    }
    stopifnot(all(rowSums(Y) == depth))
    dat.lst[[l]] <- Y
    X.lst[[l]] <- matrix(X, ncol = 1L)
    keep.lst[[l]] <- keep
    depth.lst[[l]] <- depth
    cluster.lst[[l]] <- cluster
  }
  list(OTU = dat.lst, X = X.lst, keep = keep.lst, depth = depth.lst,
       cluster = cluster.lst, correlated = correlated,
       target.taxon = taxa.names[target.idx],
       baseline.direction = baseline.direction, directions = directions)
}


# Fit each count/zero null model once, then reuse the score summaries for FE/RE.
evaluate_meta_methods <- function(simulated, sim.cnt) {
  started <- proc.time()[["elapsed"]]
  p <- statistic <- setNames(rep(NA_real_, length(meta.methods)), meta.methods)
  asymptotic <- setNames(rep(NA_real_, length(meta.methods)), meta.methods)
  reason <- warnings <- setNames(rep("", length(meta.methods)), meta.methods)
  completed <- setNames(integer(length(meta.methods)), meta.methods)
  null.seconds <- boot.seconds <- setNames(rep(0, length(meta.methods)), meta.methods)
  generation.seconds <- preprocessing.seconds <- 0
  n.samples <- n.subjects <- NA_integer_
  generation.error <- tryCatch({
    counts <- setNames(simulated$OTU, study.names)
    X <- setNames(simulated$X, study.names)
    for (s in study.names) {
      index <- match(s, study.names)
      rownames(counts[[s]]) <- rownames(X[[s]]) <-
        paste0(s, "_sample", simulated$keep[[index]])
      colnames(X[[s]]) <- "disease"
    }
    cluster.id <- if (isTRUE(simulated$correlated))
      setNames(lapply(simulated$cluster, factor), study.names) else NULL
    t0 <- proc.time()[["elapsed"]]
    d <- MetaQCAT:::.prepare_inputs(counts, X, 1L, cluster.id = cluster.id)
    ids <- MetaQCAT:::.complete_cluster_ids(d$counts, d$cluster.id)
    weights <- MetaQCAT:::.make_cluster_weights(ids, n.boot, seed + sim.cnt)
    design <- lapply(d$X, function(x) cbind(`(Intercept)` = 1, x))
    rarefied <- MetaQCAT:::.with_seed(as.integer((as.double(seed + sim.cnt) + 1) %%
      .Machine$integer.max), lapply(d$counts, MetaQCAT:::.rarefy_counts))
    zeros <- lapply(rarefied, function(y) (y == 0) * 1)
    preprocessing.seconds <- proc.time()[["elapsed"]] - t0
    n.samples <- sum(vapply(d$counts, nrow, integer(1)))
    n.subjects <- sum(vapply(ids, nlevels, integer(1)))
    NULL
  }, error = function(e) conditionMessage(e))

  diagnostics <- list()
  if (!is.null(generation.error)) reason[] <- generation.error else {
    for (component in c("count", "zero")) {
      component.methods <- paste0(c("FE-MetaQCAT", "RE-MetaQCAT"),
                                  if (component == "zero") "*" else "")
      component.warnings <- character()
      t0 <- proc.time()[["elapsed"]]
      summary <- tryCatch(withCallingHandlers({
        if (component == "count") MetaQCAT:::.Score_test_QCAT_stat(
          d$counts, design, 2L, prev.filter = min.prev, cluster.list = ids,
          meta.method = "FE-MetaQCAT") else MetaQCAT:::.Score_test_zero_stat(
          zeros, design, 2L, prev.filter = min.prev, cluster.list = ids,
          meta.method = "FE-MetaQCAT*")
      }, warning = function(w) {
        component.warnings <<- c(component.warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = identity)
      null.seconds[component.methods] <- proc.time()[["elapsed"]] - t0
      warnings[component.methods] <- paste(unique(component.warnings), collapse = " | ")
      if (inherits(summary, "error")) {
        reason[component.methods] <- conditionMessage(summary)
        next
      }
      diagnostics[[component]] <- summary$score.diagnostics
      asymptotic[component.methods[1L]] <- summary$score.pvalue
      for (method in component.methods) {
        t0 <- proc.time()[["elapsed"]]
        boot <- tryCatch(MetaQCAT:::score_bootstrap_meta(summary,
          MetaQCAT:::.weights_for_summary(summary, weights), summary$col.index.list,
          if (component == "count") summary$n.par.interest.beta else
            summary$n.par.interest.alpha, 1L, method, component), error = identity)
        boot.seconds[method] <- proc.time()[["elapsed"]] - t0
        if (inherits(boot, "error")) reason[method] <- conditionMessage(boot) else {
          value <- tryCatch(MetaQCAT:::.bootstrap_p(boot, n.boot), error = identity)
          if (inherits(value, "error")) reason[method] <- conditionMessage(value) else {
            p[method] <- value
            statistic[method] <- boot$observed.statistic
            completed[method] <- length(boot$bootstrap.statistic)
          }
        }
      }
    }
    for (model in c("FE", "RE")) {
      components <- paste0(model, c("-MetaQCAT", "-MetaQCAT*"))
      method <- paste0(model, "-MetaQCAT-O")
      null.seconds[method] <- sum(null.seconds[components])
      boot.seconds[method] <- sum(boot.seconds[components])
      warnings[method] <- paste(unique(warnings[components][nzchar(warnings[components])]), collapse = " | ")
      if (all(is.finite(p[components]))) {
        p[method] <- MetaQCAT:::ACAT(p[components])
        completed[method] <- min(completed[components])
      } else reason[method] <- paste(unique(reason[components][nzchar(reason[components])]), collapse = " | ")
      if (model == "FE" && all(is.finite(asymptotic[components])))
        asymptotic[method] <- MetaQCAT:::ACAT(asymptotic[components])
    }
  }
  valid <- is.finite(p) & p >= 0 & p <= 1
  reason[!valid & !nzchar(reason)] <- "No finite valid bootstrap p-value returned."
  list(rows = data.frame(replicate = sim.cnt, method = meta.methods,
    bootstrap_p = unname(p), asymptotic_p = unname(asymptotic),
    statistic = unname(statistic), failed = !valid,
    error = unname(reason), warning = unname(warnings),
    n_boot_requested = n.boot, n_boot_completed = unname(completed), n_perm = 0L,
    null_fit_seconds = unname(null.seconds), bootstrap_seconds = unname(boot.seconds),
    method_runtime_seconds = unname(null.seconds + boot.seconds),
    generation_seconds = generation.seconds, preprocessing_seconds = preprocessing.seconds,
    paired_replicate_runtime_seconds = proc.time()[["elapsed"]] - started,
    n_samples = n.samples, n_subjects = n.subjects,
    data_seed = seed + sim.cnt, bootstrap_seed = seed + sim.cnt,
    target_taxon = if (is.null(simulated)) NA_character_ else simulated$target.taxon,
    directions = if (is.null(simulated)) NA_character_ else paste(simulated$directions, collapse = ";"),
    stringsAsFactors = FALSE), diagnostics = diagnostics)
}

# Comparator rarefaction and Fisher combination used in the paper scripts.
Rarefy <- function(otu.tab, depth = min(rowSums(otu.tab)), seed = 1234L) {
  otu.tab <- as.matrix(otu.tab)
  keep <- rowSums(otu.tab) >= depth
  otu.tab <- otu.tab[keep, , drop = FALSE]
  set.seed(seed)
  out <- t(apply(otu.tab, 1L, function(x) {
    draw <- sample(rep(seq_along(x), x), depth)
    tabulate(draw, nbins = length(x))
  }))
  dimnames(out) <- dimnames(otu.tab)
  out
}

F.test <- function(x) {
  x <- x[is.finite(x) & x >= 0 & x <= 1]
  if (!length(x)) return(NA_real_)
  if (any(x == 0)) return(0)
  stats::pchisq(-2 * sum(log(x)), df = 2 * length(x), lower.tail = FALSE)
}

# DM preserves the comparator's exceedances / successful-permutations rule.
# Unlike MetaQCAT, this comparator does not use the plus-one correction.
dm1.resample <- function(Y, case, n.perm) {
  observed <- HMP::Xmcupo.sevsample(list(
    grp1 = Y[case == 0, , drop = FALSE],
    grp2 = Y[case == 1, , drop = FALSE]))
  statistic <- as.numeric(observed[["Xmcupo statistics"]])
  if (length(statistic) != 1L || !is.finite(statistic))
    stop("DM returned an invalid observed statistic.")
  completed <- exceed <- 0L
  for (b in seq_len(n.perm)) {
    label <- sample(case)
    value <- tryCatch(HMP::Xmcupo.sevsample(list(
      grp1 = Y[label == 0, , drop = FALSE],
      grp2 = Y[label == 1, , drop = FALSE]))[["Xmcupo statistics"]],
      error = function(e) NA_real_)
    if (length(value) == 1L && is.finite(value)) {
      completed <- completed + 1L
      exceed <- exceed + (value >= statistic)
    }
  }
  list(asymptotic = as.numeric(observed[["p value"]]),
       resampling = if (completed > 0L) exceed / completed else NA_real_,
       completed = completed)
}

# No data generation occurs in this function. P uses single-study tests with
# Fisher combination; CD adjusts pooled data with MMUPHin; RD pools raw data.
evaluate_comparators <- function(simulated, sim.cnt) {
  output <- vector("list", 12L)
  at <- 0L
  for (strategy in c("P", "CD", "RD")) {
    set.seed(seed + sim.cnt)
    preparation.started <- proc.time()[["elapsed"]]
    tables <- simulated$OTU
    designs <- simulated$X
    preparation.error <- ""
    if (strategy != "P") {
      tables <- list(do.call(rbind, tables))
      designs <- list(matrix(unlist(designs, use.names = FALSE), ncol = 1L))
    }
    if (strategy == "CD") {
      corrected <- tryCatch({
        metadata <- data.frame(
          Country = factor(rep(study.names, lengths(simulated$keep)),
                           levels = study.names),
          Group = factor(ifelse(designs[[1L]][, 1L] == 1, "case", "control")))
        rownames(metadata) <- paste0("sample", seq_len(nrow(metadata)))
        abundance <- t(tables[[1L]])
        colnames(abundance) <- rownames(metadata)
        MMUPHin::adjust_batch(
          feature_abd = abundance, batch = "Country", covariates = "Group",
          data = metadata, control = list(verbose = FALSE))$feature_abd_adj
      }, error = identity)
      if (inherits(corrected, "error")) {
        preparation.error <- conditionMessage(corrected)
      } else tables[[1L]] <- t(corrected)
    }
    rarefied <- if (!nzchar(preparation.error)) tryCatch(
      lapply(tables, Rarefy), error = identity) else NULL
    if (inherits(rarefied, "error"))
      preparation.error <- conditionMessage(rarefied)
    preparation.seconds <- proc.time()[["elapsed"]] - preparation.started

    for (test in c("QCAT1", "QCAT2", "DM", "ZIGDM")) {
      at <- at + 1L
      budget <- if (strategy == "RD" && test == "DM") n.perm.rd.dm else n.perm
      p <- asymptotic <- rep(NA_real_, length(tables))
      completed <- rep(NA_integer_, length(tables))
      errors <- warnings <- character()
      started <- proc.time()[["elapsed"]]
      if (!nzchar(preparation.error)) for (s in seq_along(tables)) {
        value <- tryCatch(withCallingHandlers({
          Y <- tables[[s]]
          X <- designs[[s]]
          if (test == "DM") {
            dm1.resample(Y, X[, 1L], budget)
          } else {
            fit <- switch(test,
              QCAT1 = miLineage::QCAT(Y, X, 1L, Tax = NULL,
                                      min.depth = 0, n.perm = budget),
              QCAT2 = miLineage::QCAT_GEE(rarefied[[s]], X, 1L, X, 1L,
                                          Tax = NULL, min.depth = 0,
                                          n.perm = budget),
              ZIGDM = miLineage::ZIGDM(Y, X, NULL, NULL, Tax = NULL,
                                        X.index = 1L, n.perm = budget))
            list(asymptotic = as.numeric(fit$pval[1L]),
                 resampling = as.numeric(fit$pval[2L]),
                 completed = NA_integer_)
          }
        }, warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }, message = function(m) invokeRestart("muffleMessage")),
        error = identity)
        if (inherits(value, "error")) {
          errors <- c(errors, paste0("study ", s, ": ", conditionMessage(value)))
        } else {
          p[s] <- value$resampling
          asymptotic[s] <- value$asymptotic
          completed[s] <- value$completed
        }
      }
      valid.study <- is.finite(p) & p >= 0 & p <= 1
      combined <- if (strategy == "P") F.test(p) else p[1L]
      combined.asy <- if (strategy == "P") F.test(asymptotic) else asymptotic[1L]
      failed <- !is.finite(combined) || combined < 0 || combined > 1
      reason <- paste(c(preparation.error, unique(errors))[
        nzchar(c(preparation.error, unique(errors)))], collapse = " | ")
      if (failed && !nzchar(reason)) reason <- "No valid resampling p-value."
      output[[at]] <- data.frame(
        method = paste0(strategy, "-", test, " (perm)"),
        resampling_p = combined, asymptotic_p = combined.asy, failed = failed,
        error = reason, warning = paste(unique(warnings), collapse = " | "),
        n_boot_requested = 0L, n_boot_completed = 0L,
        n_perm_configured = budget,
        n_perm_completed = if (test == "DM" && any(!is.na(completed)))
          sum(completed, na.rm = TRUE) else NA_integer_,
        n_perm_tests = length(tables), n_studies_used = sum(valid.study),
        partial_study_failure = any(!valid.study),
        method_runtime_seconds = proc.time()[["elapsed"]] - started +
          preparation.seconds,
        null_fit_seconds = NA_real_, bootstrap_seconds = NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, output)
}

# Native cluster-QCAT calls for correlated null and power settings.
# BTW internally concatenates clusters, so sort Y/X/ID jointly before calling.
# All-zero composition rows are retained for the observed-zero component.
run_cluster_qcat <- function(y, x, id, two.part = FALSE) {
  warnings <- character()
  fit <- tryCatch(withCallingHandlers({
    y <- as.matrix(y)
    x <- as.matrix(x)
    id <- as.character(id)
    if (!nrow(y)) stop("No samples.")
    if (nrow(y) != nrow(x) || nrow(y) != length(id))
      stop("Count rows, covariate rows and cluster IDs must align.")
    order.id <- order(id, rownames(y), method = "radix")
    y <- y[order.id, , drop = FALSE]
    x <- x[order.id, , drop = FALSE]
    id <- factor(id[order.id])
    if (anyNA(x) || anyNA(id) || length(unique(x[, 1L])) != 2L)
      stop("Disease must be observed and variable after sample screening.")
    if (any(vapply(split(x[, 1L], id), function(z)
      length(unique(z)) != 1L, logical(1))))
      stop("BTW permutation requires disease to be constant within each cluster.")
    if (sum(colSums(y) > 0) < 2L) stop("Fewer than two observed taxa.")
    if (two.part) {
      miLineage::QCAT_GEE.Cluster(
        ID = id, OTU = y, X = x, X.index = 1L, Z = x, Z.index = 1L,
        Tax = NULL, min.depth = 0, perm.type = "BTW",
        n.perm = n.perm, test = "chisq")
    } else {
      miLineage::QCAT.Cluster(
        ID = id, OTU = y, X = x, X.index = 1L,
        Tax = NULL, min.depth = 0, perm.type = "BTW",
        n.perm = n.perm, test = "chisq")
    }
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")), error = identity)
  p <- c(Asymptotic = NA_real_, Resampling = NA_real_)
  reason <- ""
  if (inherits(fit, "error")) reason <- conditionMessage(fit) else {
    value <- tryCatch(unlist(fit$pval), error = function(e) NULL)
    if (!is.numeric(value) || is.null(names(value))) {
      reason <- "Native QCAT returned no named numeric p-values."
    } else {
      p[] <- value[match(names(p), names(value))]
      if (any(!is.finite(p))) reason <- "Native QCAT returned nonfinite p-values."
    }
  }
  bad <- !is.finite(p) | p < 0 | p > 1
  p[bad] <- NA_real_
  list(p = p, error = reason,
       warning = paste(unique(warnings), collapse = " | "))
}

# Same result schema as evaluate_comparators(), but only six cluster-aware
# QCAT tests. No generation occurs here; independent tests stay above.
evaluate_cluster_comparators <- function(simulated, sim.cnt) {
  if (!isTRUE(simulated$correlated))
    stop("Cluster comparators require the correlated-sample generator.")
  counts <- setNames(simulated$OTU, study.names)
  designs <- setNames(simulated$X, study.names)
  ids <- setNames(lapply(simulated$cluster, factor), study.names)
  for (s in seq_along(counts)) {
    rownames(counts[[s]]) <- rownames(designs[[s]]) <-
      paste0(study.names[s], "_sample", simulated$keep[[s]])
    colnames(designs[[s]]) <- "disease"
    if (length(ids[[s]]) != nrow(counts[[s]]))
      stop("Generated cluster IDs do not align with count rows.")
  }
  pooled <- do.call(rbind, unname(counts))
  pooled.X <- do.call(rbind, unname(designs))
  # Prefix cluster IDs to prevent distinct studies sharing a pooled subject.
  pooled.ids <- factor(unlist(Map(function(study, id)
    paste(study, as.character(id), sep = "::"), names(ids), ids),
    use.names = FALSE))
  metadata <- data.frame(
    Country = factor(rep(study.names, vapply(counts, nrow, integer(1)))),
    Group = factor(ifelse(pooled.X[, 1L] == 1, "case", "control")),
    row.names = rownames(pooled))
  output <- vector("list", 6L)
  at <- 0L
  for (strategy in c("P", "CD", "RD")) {
    set.seed(seed + sim.cnt)
    preparation.started <- proc.time()[["elapsed"]]
    tables <- if (strategy == "P") counts else list(pooled)
    X <- if (strategy == "P") designs else list(pooled.X)
    clusters <- if (strategy == "P") ids else list(pooled.ids)
    preparation.error <- ""
    if (strategy == "CD") {
      corrected <- tryCatch(MMUPHin::adjust_batch(
        feature_abd = t(pooled), batch = "Country", covariates = "Group",
        data = metadata,
        control = list(verbose = FALSE, diagnostic_plot = NULL))$feature_abd_adj,
        error = identity)
      if (inherits(corrected, "error"))
        preparation.error <- conditionMessage(corrected) else
        tables[[1L]] <- t(corrected)
    }
    preparation.seconds <- proc.time()[["elapsed"]] - preparation.started

    # Raw counts for QCAT.Cluster; rarefy each study (P) or pooled table (CD/RD)
    # once before QCAT_GEE.Cluster. Keep the generated rows and IDs aligned.
    rarefaction.started <- proc.time()[["elapsed"]]
    set.seed(seed + sim.cnt)
    rarefied <- if (!nzchar(preparation.error))
      tryCatch(lapply(tables, MetaQCAT:::.rarefy_counts), error = identity) else NULL
    rarefaction.error <- if (inherits(rarefied, "error"))
      conditionMessage(rarefied) else ""
    rarefaction.seconds <- proc.time()[["elapsed"]] - rarefaction.started
    for (test in c("QCAT.Cluster", "QCAT_GEE.Cluster")) {
      at <- at + 1L
      two.part <- test == "QCAT_GEE.Cluster"
      p.study <- matrix(NA_real_, length(tables), 2L,
        dimnames = list(names(tables), c("Asymptotic", "Resampling")))
      errors <- warnings <- character()
      started <- proc.time()[["elapsed"]]
      if (nzchar(preparation.error)) errors <- preparation.error else if (
        two.part && nzchar(rarefaction.error)) errors <- rarefaction.error else {
        for (s in seq_along(tables)) {
          value <- run_cluster_qcat(
            if (two.part) rarefied[[s]] else tables[[s]],
            X[[s]], clusters[[s]], two.part = two.part)
          p.study[s, ] <- value$p
          if (nzchar(value$error))
            errors <- c(errors, paste0("study ", s, ": ", value$error))
          if (nzchar(value$warning)) warnings <- c(warnings, value$warning)
        }
      }
      valid.study <- is.finite(p.study[, "Resampling"])
      combined <- if (strategy == "P") F.test(p.study[, "Resampling"]) else
        p.study[1L, "Resampling"]
      combined.asy <- if (strategy == "P") F.test(p.study[, "Asymptotic"]) else
        p.study[1L, "Asymptotic"]
      failed <- !is.finite(combined) || combined < 0 || combined > 1
      reason <- paste(unique(errors), collapse = " | ")
      if (failed && !nzchar(reason)) reason <- "No valid resampling p-value."
      output[[at]] <- data.frame(
        method = paste0(strategy, "-", test, " (perm)"),
        resampling_p = combined, asymptotic_p = combined.asy, failed = failed,
        error = reason, warning = paste(unique(warnings), collapse = " | "),
        n_boot_requested = 0L, n_boot_completed = 0L,
        n_perm_configured = n.perm,
        # Native adaptive completed permutation counts are not returned.
        n_perm_completed = NA_integer_, n_perm_tests = length(tables),
        n_studies_used = sum(valid.study),
        partial_study_failure = strategy == "P" && any(valid.study) && any(!valid.study),
        method_runtime_seconds = proc.time()[["elapsed"]] - started +
          preparation.seconds + if (two.part) rarefaction.seconds else 0,
        null_fit_seconds = NA_real_, bootstrap_seconds = NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  # P calls permute within studies; native CD/RD BTW calls permute across the
  # pooled clusters. Native functions do not accept an extra study block.
  do.call(rbind, output)
}

# The only generation call in the complete all-method replicate workflow.
run_simulation_replicate <- function(sim.cnt, signal.size, signal.method,
                                     type, direct, correlated = FALSE) {
  started <- proc.time()[["elapsed"]]
  cluster.mode <- isTRUE(correlated)
  methods <- c(meta.methods, if (cluster.mode)
    cluster.baseline.methods else baseline.methods)
  rows <- data.frame(
    method = methods, resampling_p = NA_real_, asymptotic_p = NA_real_,
    failed = TRUE, error = "", warning = "", n_boot_requested = 0L,
    n_boot_completed = 0L, n_perm_configured = 0L,
    n_perm_completed = NA_integer_, n_perm_tests = 0L,
    n_studies_used = 0L, partial_study_failure = FALSE,
    method_runtime_seconds = 0, null_fit_seconds = NA_real_,
    bootstrap_seconds = NA_real_, stringsAsFactors = FALSE)
  rows$n_boot_requested[match(meta.methods, rows$method)] <- n.boot
  baseline.index <- !rows$method %in% meta.methods
  rows$n_perm_configured[baseline.index] <- n.perm
  if (!cluster.mode)
    rows$n_perm_configured[rows$method == "RD-DM (perm)"] <- n.perm.rd.dm
  simulated <- tryCatch(simulate_crc_data(
    sim.cnt, mod_gdm_lst, depth.pool, taxa.names, lst.num,
    signal.size, signal.method, type, direct, seed, correlated),
    error = identity)
  generation.seconds <- proc.time()[["elapsed"]] - started
  diagnostics <- list()
  if (inherits(simulated, "error")) {
    rows$error <- paste("Data generation:", conditionMessage(simulated))
  } else {
    fit <- evaluate_meta_methods(simulated, sim.cnt)
    diagnostics <- fit$diagnostics
    result <- fit$rows
    result$resampling_p <- result$bootstrap_p
    result$n_perm_configured <- 0L
    result$n_perm_completed <- 0L
    result$n_perm_tests <- 0L
    result$n_studies_used <- ifelse(result$failed, 0L, length(study.names))
    common <- intersect(names(rows), names(result))
    index <- match(result$method, rows$method)
    if (anyNA(index) || anyDuplicated(index) ||
        !setequal(result$method, meta.methods))
      stop("MetaQCAT method rows do not match the configured simulation methods.")
    rows[index, common] <- result[common]
    comparator <- if (cluster.mode)
      evaluate_cluster_comparators(simulated, sim.cnt) else
      evaluate_comparators(simulated, sim.cnt)
    index <- match(comparator$method, rows$method)
    expected <- if (cluster.mode) cluster.baseline.methods else baseline.methods
    if (anyNA(index) || anyDuplicated(index) ||
        !setequal(comparator$method, expected))
      stop("Comparator rows do not match the configured simulation methods.")
    rows[index, names(comparator)] <- comparator
  }
  rows$replicate <- sim.cnt
  rows$calibration <- ifelse(rows$method %in% meta.methods,
    "wild-score bootstrap", if (cluster.mode)
      "between-cluster permutation" else "permutation")
  rows$generation_seconds <- generation.seconds
  rows$replicate_runtime_seconds <- proc.time()[["elapsed"]] - started
  rows$n_samples <- if (inherits(simulated, "error")) NA_integer_ else
    sum(vapply(simulated$OTU, nrow, integer(1)))
  rows$n_subjects <- if (inherits(simulated, "error")) NA_integer_ else
    sum(vapply(simulated$cluster, function(x) length(unique(x)), integer(1)))
  rows$data_seed <- seed + sim.cnt
  rows$bootstrap_seed <- ifelse(rows$method %in% meta.methods, seed + sim.cnt, NA)
  rows$target_taxon <- if (inherits(simulated, "error")) NA_character_ else
    simulated$target.taxon
  rows$directions <- if (inherits(simulated, "error")) NA_character_ else
    paste(simulated$directions, collapse = ";")
  list(rows = rows, diagnostics = diagnostics)
}
