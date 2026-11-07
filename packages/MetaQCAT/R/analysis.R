# Rarefy once per study before any lineage aggregation, as in MetaPALM.
.rarefy_counts <- function(counts) {
  depth <- min(rowSums(counts))
  draw <- function(x) {
    out <- numeric(length(x)); remaining <- depth; population <- sum(x)
    if (length(x) > 1L) {
      for (j in seq_len(length(x) - 1L)) {
        out[j] <- rhyper(1L, x[j], population - x[j], remaining)
        remaining <- remaining - out[j]; population <- population - x[j]
      }
    }
    out[length(x)] <- remaining
    out
  }
  out <- t(vapply(seq_len(nrow(counts)), function(i) draw(counts[i, ]),
                  numeric(ncol(counts))))
  dimnames(out) <- dimnames(counts)
  out
}

.with_seed <- function(seed, expr) {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old.kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old.kind))
    if (existed) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  force(expr)
}

.complete_cluster_ids <- function(counts, cluster.id) {
  if (is.null(cluster.id)) cluster.id <- vector("list", length(counts))
  if (is.null(names(counts))) names(counts) <- paste0("Study", seq_along(counts))
  names(cluster.id) <- names(counts)
  for (i in seq_along(cluster.id)) {
    id <- if (is.null(cluster.id[[i]])) {
      if (is.null(rownames(counts[[i]]))) paste0("Sample", seq_len(nrow(counts[[i]])))
      else rownames(counts[[i]])
    } else as.character(cluster.id[[i]])
    if (anyNA(id) || any(!nzchar(id))) stop("Subject IDs cannot be missing or empty.")
    cluster.id[[i]] <- factor(id, levels = sort(unique(id), method = "radix"))
  }
  cluster.id
}

.make_cluster_weights <- function(cluster.id, n.boot, seed) {
  .check_scalar(n.boot, "n.boot", 1, .Machine$integer.max, integer = TRUE)
  if (is.null(names(cluster.id))) names(cluster.id) <- paste0("Study", seq_along(cluster.id))
  if (anyNA(names(cluster.id)) || any(!nzchar(names(cluster.id))) ||
      anyDuplicated(names(cluster.id))) stop("Unique study keys are required.")
  # Generate in canonical key order, then return in the caller's study order.
  # Filtering later selects rows by ID rather than shifting the weight plan.
  keys <- sort(names(cluster.id), method = "radix")
  out <- .with_seed(seed, lapply(cluster.id[keys], function(id) {
    lev <- sort(unique(as.character(id)), method = "radix")
    if (anyNA(lev) || any(!nzchar(lev))) stop("Subject IDs cannot be missing or empty.")
    matrix(sample(c(-1, 1), length(lev) * n.boot, replace = TRUE),
           nrow = length(lev), ncol = n.boot, dimnames = list(lev, NULL))
  }))
  out[names(cluster.id)]
}

.weights_for_summary <- function(summary, weights) {
  names.included <- names(summary$cluster.list)
  if (is.null(names.included) || any(!names.included %in% names(weights)))
    stop("Study keys were lost during component filtering.")
  lapply(names.included, function(name) {
    id <- summary$cluster.list[[name]]
    keys <- levels(droplevels(id))
    w <- weights[[name]]
    if (!all(keys %in% rownames(w)))
      stop("Cluster keys were lost during component filtering.")
    w[keys, , drop = FALSE]
  })
}

.bootstrap_metadata <- function(weights, seed, completed) {
  list(calibration = "null-imposed wild-score bootstrap",
       multipliers = "subject-level Rademacher", centered = FALSE,
       covariance = "fixed unscaled cluster sandwich",
       seed = as.integer(seed), requested = as.integer(ncol(weights[[1L]])),
       completed = as.integer(completed),
       rng = c("Mersenne-Twister", "Inversion", "Rejection"),
       subject.ids = lapply(weights, rownames),
       n.clusters = vapply(weights, nrow, integer(1)))
}

.attach_score_keys <- function(res, ids, maps, taxa, targets) {
  studies <- names(ids)
  if (is.null(studies)) studies <- paste0("Study", seq_along(ids))
  coef.ids <- as.vector(t(outer(taxa, targets, paste, sep = ":")))
  fields <- intersect(c("psi.lst", "cov.lst", "bread.lst", "sensitivity.lst",
                        "cluster.score.lst", "score.diagnostics",
                        "kernel.lst", "covariance.diagnostics"), names(res))
  for (field in fields) names(res[[field]]) <- studies
  for (i in seq_along(res$score.diagnostics)) {
    if (!is.finite(res$score.diagnostics[[i]]$nuisance.score.relative) ||
        res$score.diagnostics[[i]]$nuisance.score.relative > 1e-8)
      stop("Null nuisance equations were not solved accurately in study ", studies[i],
           "; relative nuisance-score residual exceeds 1e-8.")
  }
  subject.ids <- lapply(ids, function(id) levels(droplevels(id)))
  if (is.null(res$psi.lst)) res$psi.lst <- lapply(seq_along(ids), function(i)
    t(res$bread.lst[[i]] %*% res$cluster.score.lst[[i]]))
  names(res$psi.lst) <- names(subject.ids) <- names(maps) <- studies
  for (i in seq_along(ids)) {
    local.ids <- coef.ids[maps[[i]]]
    dimnames(res$psi.lst[[i]]) <- list(subject.ids[[i]], local.ids)
    dimnames(res$cov.lst[[i]]) <- list(local.ids, local.ids)
  }
  res$subject.ids <- subject.ids
  res$coefficient.ids <- coef.ids
  res$alignment <- maps
  res$study.ids <- studies
  res
}

.component_record <- function(summary) {
  summary[intersect(c("psi.lst", "cov.lst", "sensitivity.lst", "subject.ids",
    "coefficient.ids", "alignment", "study.ids", "score.diagnostics",
    "kernel.lst", "covariance.diagnostics", "reference.taxon", "re.diagnostics"), names(summary))]
}

.bootstrap_p <- function(boot, n.boot) {
  if (length(boot$bootstrap.statistic) != n.boot ||
      any(!is.finite(c(boot$observed.statistic, boot$bootstrap.statistic))))
    stop("Every requested bootstrap statistic must be finite; no draws may be discarded.")
  (1 + sum(boot$bootstrap.statistic >= boot$observed.statistic)) / (n.boot + 1)
}

.run_meta_test <- function(counts, zeros, X, X.index, cluster.id, weights, method,
                            min.prev, n.boot, save.coef, seed,
                            reference.taxon = NULL) {
  args <- list(X.list = X, X.par.index = X.index, prev.filter = min.prev[1L],
               taxa.names = NULL, cluster.list = cluster.id, weights = weights,
               n.boot = n.boot,
               meta.method = sub("-O$", "", method), save.desc = save.coef,
               seed = seed)
  composition <- function() do.call(.Score_test_QCAT_meta,
                                    c(list(Y.list = counts,
                                           reference.taxon = reference.taxon), args))
  zero <- function() {
    do.call(.Score_test_zero_meta, list(
      Y.bin.list = zeros, Z.list = X, Z.par.index = X.index,
      prev.filter = min.prev[length(min.prev)], taxa.names = NULL,
      cluster.list = cluster.id, weights = weights, n.boot = n.boot,
      meta.method = if (grepl("-O$", method)) paste0(sub("-O$", "", method), "*") else method,
      save.desc = save.coef, seed = seed))
  }
  tryCatch({
    if (grepl("\\*$", method)) {
      res <- zero()
    } else if (grepl("-O$", method)) {
      res <- composition(); z <- zero()
      if (!identical(names(res$pvalue.comp), names(z$pvalue.pa)))
        stop("Component calibrations do not align.")
      combined <- c(Boot = ACAT(c(res$pvalue.comp["Boot"], z$pvalue.pa["Boot"])))
      if ("Asy" %in% names(res$pvalue.comp))
        combined <- c(Asy = ACAT(c(res$pvalue.comp["Asy"], z$pvalue.pa["Asy"])),
                      combined)
      eligibility <- c(res$eligibility, z$eligibility)
      res <- c(res, z[setdiff(names(z), c("bootstrap", "eligibility", "status"))],
               list(pvalue.comb = combined,
                    bootstrap.omnibus = list(
                      calibration = "Cauchy combination of marginal bootstrap p-values",
                      weights = c(count = .5, zero = .5))))
      res$eligibility <- eligibility
    } else {
      res <- composition()
    }
    p <- .extract_pvalues(res, method)
    if (!length(p) || any(!is.finite(p)) || any(p < 0 | p > 1))
      stop("The backend returned invalid p-values.")
    res
  }, error = function(e) {
    calibration <- if (grepl("^RE-", method)) "Boot" else c("Asy", "Boot")
    field <- if (grepl("-O$", method)) "pvalue.comb" else
      if (grepl("\\*$", method)) "pvalue.pa" else "pvalue.comp"
    res <- list(error = conditionMessage(e))
    res$status <- "nonestimable"
    res[[field]] <- setNames(rep(NA_real_, length(calibration)), calibration)
    res$bootstrap <- .bootstrap_metadata(weights, seed, completed = 0L)
    res
  })
}

.extract_pvalues <- function(result, method) {
  field <- if (grepl("-O$", method)) "pvalue.comb" else
    if (grepl("\\*$", method)) "pvalue.pa" else "pvalue.comp"
  result[[field]]
}

.scan_lineages <- function(counts, zeros, X, X.index, cluster.id, weights, taxonomy,
                            method, min.prev, n.boot, fdr.alpha,
                            save.coef, n.threads, seed) {
  tasks <- list()
  for (rank in rev(seq_len(ncol(taxonomy) - 1L))) {
    for (parent in unique(taxonomy[, rank])) {
      rows <- which(taxonomy[, rank] == parent)
      children <- taxonomy[rows, rank + 1L]
      if (length(unique(children)) > 1L)
        tasks[[length(tasks) + 1L]] <- list(parent = parent, rank = rank,
                                          rows = rows, children = children)
    }
  }
  keys <- vapply(tasks, function(t) t$parent, character(1))
  duplicate <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
  keys[duplicate] <- vapply(tasks[duplicate], function(t)
    paste0("Rank", t$rank, ":", t$parent), character(1))
  keys <- make.unique(keys)
  aggregate <- function(y, task) {
    if (is.null(y)) return(NULL)
    lapply(y, function(m) {
      a <- t(rowsum(t(m[, task$rows, drop = FALSE]), task$children, reorder = FALSE))
      rownames(a) <- rownames(m)
      a
    })
  }
  work <- function(i) {
    task <- tasks[[i]]
    y <- aggregate(counts, task)
    z <- aggregate(zeros, task)
    if (!is.null(z)) z <- lapply(z, function(m) (m == 0) * 1)
    .run_meta_test(counts = y, zeros = z, X = X, X.index = X.index,
                   cluster.id = cluster.id, weights = weights,
                   method = method, min.prev = min.prev, n.boot = n.boot,
                   save.coef = save.coef, seed = seed)
  }
  if (.Platform$OS.type == "windows" && n.threads > 1L) {
    warning("Fork-based parallelism is unavailable on Windows; using one worker.",
            call. = FALSE)
    n.threads <- 1L
  }
  results <- if (n.threads == 1L || length(tasks) < 2L) lapply(seq_along(tasks), work) else
    parallel::mclapply(seq_along(tasks), work, mc.set.seed = FALSE,
                       mc.cores = min(n.threads, parallelly::availableCores()))
  names(results) <- keys
  calibration <- if (grepl("^RE-", method)) "Boot" else c("Asy", "Boot")
  p <- matrix(NA_real_, length(calibration), length(tasks),
              dimnames = list(paste0(ifelse(calibration == "Asy", "Asymptotic-", "Bootstrap-"), method), keys))
  for (i in seq_along(results)) p[, i] <- .extract_pvalues(results[[i]], method)[calibration]
  failed <- vapply(results, function(r) if (is.null(r$error)) NA_character_ else r$error,
                   character(1))
  failed <- failed[!is.na(failed)]
  if (length(failed)) warning(length(failed),
    " lineage test(s) failed; NA columns retained. See failed.lineage for reasons.", call. = FALSE)
  q <- p
  for (i in seq_len(nrow(q))) q[i, ] <- p.adjust(p[i, ], method = "BH")
  final.q <- q[nrow(q), ]
  list(lineage.pval = p, lineage.qval = q, lineage.coef = results,
       sig.lineage = keys[!is.na(final.q) & final.q <= fdr.alpha],
       failed.lineage = failed,
       bootstrap = c(.bootstrap_metadata(weights, seed,
         if (length(results) && !length(failed)) n.boot else 0L),
         list(completed.by.lineage = vapply(results, function(r)
           r$bootstrap$completed, integer(1)))))
}
