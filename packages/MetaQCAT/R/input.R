.check_arg_names <- function(call, allowed) {
  supplied <- names(as.list(call))[-1L]
  bad <- supplied[nzchar(supplied) & !supplied %in% allowed]
  if (length(bad)) stop("unused argument name(s): ", paste(bad, collapse = ", "),
                         ". Use exact release argument names.", call. = FALSE)
}

.check_scalar <- function(x, name, lower = 0, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    stop(name, " must be a finite ", if (integer) "integer" else "number",
         " between ", lower, " and ", upper, ".", call. = FALSE)
  }
}

.check_options <- function(method, min.prev, n.boot, save.coef, seed,
                           omnibus = TRUE) {
  methods <- c("FE-MetaQCAT", "RE-MetaQCAT", "FE-MetaQCAT*", "RE-MetaQCAT*")
  if (omnibus) methods <- c(methods, "FE-MetaQCAT-O", "RE-MetaQCAT-O")
  if (!is.character(method) || length(method) != 1L || !method %in% methods)
    stop("Unknown method.", call. = FALSE)
  if (!is.numeric(min.prev) || !length(min.prev) %in% c(1L, 2L) ||
      any(!is.finite(min.prev)) || any(min.prev < 0 | min.prev > 1) ||
      (length(min.prev) == 2L && !grepl("-O$", method)))
    stop("min.prev must be in [0, 1]; two values are allowed only for omnibus tests.",
         call. = FALSE)
  .check_scalar(n.boot, "n.boot", 1, .Machine$integer.max,
                integer = TRUE)
  if (!is.logical(save.coef) || length(save.coef) != 1L || is.na(save.coef))
    stop("save.coef must be TRUE or FALSE.", call. = FALSE)
  .check_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
}

.align_studies <- function(x, study.names, what) {
  if (!is.list(x) || length(x) != length(study.names))
    stop(what, " must be a list with one element per study.", call. = FALSE)
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x)) || !setequal(names(x), study.names))
      stop(what, " study names must match counts.", call. = FALSE)
    x <- x[study.names]
  }
  names(x) <- study.names
  x
}

.prepare_inputs <- function(counts, X, X.index, cluster.id = NULL,
                             taxonomy = NULL, binary = FALSE,
                             intercept = FALSE) {
  if (!is.list(counts) || !length(counts))
    stop("counts must be a nonempty list of matrices.", call. = FALSE)
  study.names <- names(counts)
  if (is.null(study.names)) study.names <- paste0("Study", seq_along(counts))
  if (anyNA(study.names) || any(!nzchar(study.names)) || anyDuplicated(study.names))
    stop("Study names must be nonmissing and unique.", call. = FALSE)
  names(counts) <- study.names
  X <- .align_studies(X, study.names, "X")
  if (!is.null(cluster.id)) cluster.id <- .align_studies(cluster.id, study.names,
                                                        "cluster.id")
  if (!is.numeric(X.index) || !length(X.index) || any(!is.finite(X.index)) ||
      any(X.index != floor(X.index)) || any(X.index < 1) || anyDuplicated(X.index))
    stop("X.index must contain unique positive column indices.", call. = FALSE)
  X.index <- as.integer(X.index)
  if (intercept && 1L %in% X.index)
    stop("The low-level interface cannot test the intercept (column 1).", call. = FALSE)
  for (i in seq_along(counts)) {
    y <- as.matrix(counts[[i]]); x <- as.matrix(X[[i]])
    if (!is.numeric(y) || !nrow(y) || !ncol(y) || any(!is.finite(y)) ||
        any(y < 0) || any(abs(y - round(y)) > 1e-7) ||
        (binary && any(!y %in% c(0, 1))))
      stop("Invalid ", if (binary) "zero-indicator" else "count",
           " matrix in study ", study.names[i], ".", call. = FALSE)
    if (is.null(colnames(y)) || anyNA(colnames(y)) ||
        any(!nzchar(colnames(y))) || anyDuplicated(colnames(y)))
      stop("Taxon column names must be nonmissing and unique.", call. = FALSE)
    if (!is.numeric(x) || !ncol(x) || nrow(x) != nrow(y) ||
        any(!is.finite(x)) || any(X.index > ncol(x)))
      stop("Invalid X or X.index in study ", study.names[i], ".", call. = FALSE)
    if (!is.null(rownames(y)) && !is.null(rownames(x)) &&
        !identical(rownames(y), rownames(x)))
      stop("Sample row names in counts and X do not align in study ",
           study.names[i], ".", call. = FALSE)
    if (intercept && any(x[, 1L] != 1))
      stop("Low-level X matrices must have an intercept of ones in column 1.",
           call. = FALSE)
    checked.x <- if (intercept) x else cbind(1, x)
    if (qr(checked.x)$rank < ncol(checked.x))
      stop("X is rank deficient in study ", study.names[i],
           "; remove constant or redundant covariates (do not add an intercept to MetaQCAT).",
           call. = FALSE)
    if (!is.null(cluster.id) && !is.null(cluster.id[[i]])) {
      id <- cluster.id[[i]]
      if (!is.atomic(id) || length(id) != nrow(y) || anyNA(id))
        stop("cluster.id must have one nonmissing ID per sample.", call. = FALSE)
      if (!is.null(names(id)) && !is.null(rownames(y)) &&
          !identical(names(id), rownames(y)))
        stop("Named cluster.id entries do not align with samples.", call. = FALSE)
      id <- as.character(id)
      if (any(!nzchar(id))) stop("cluster.id cannot contain empty IDs.", call. = FALSE)
      id <- factor(id, levels = sort(unique(id), method = "radix"))
      cluster.id[[i]] <- id
    }
    sample.ids <- if (!is.null(rownames(y))) rownames(y) else rownames(x)
    if (is.null(sample.ids)) sample.ids <- paste0("Sample", seq_len(nrow(y)))
    if (anyNA(sample.ids) || any(!nzchar(sample.ids)) || anyDuplicated(sample.ids))
      stop("Sample row names must be nonmissing and unique within study.", call. = FALSE)
    rownames(y) <- rownames(x) <- sample.ids
    at <- order(sample.ids, method = "radix")
    y <- y[at, , drop = FALSE]; x <- x[at, , drop = FALSE]
    if (!is.null(cluster.id) && !is.null(cluster.id[[i]]))
      cluster.id[[i]] <- cluster.id[[i]][at]
    storage.mode(y) <- storage.mode(x) <- "double"
    counts[[i]] <- y; X[[i]] <- x
  }
  target.names <- lapply(X, function(x) colnames(x)[X.index])
  if (all(vapply(target.names, length, integer(1)) == length(X.index)) &&
      any(!vapply(target.names, identical, logical(1), target.names[[1L]])))
    stop("Tested covariate names and ordering must agree across studies.", call. = FALSE)
  counts <- counts[sort(study.names, method = "radix")]
  X <- X[names(counts)]
  if (!is.null(cluster.id)) cluster.id <- cluster.id[names(counts)]
  taxa <- sort(unique(unlist(lapply(counts, colnames), use.names = FALSE)), method = "radix")
  counts <- lapply(counts, function(y) {
    if (identical(colnames(y), taxa)) return(y)
    z <- matrix(if (binary) 1 else 0, nrow(y), length(taxa),
                dimnames = list(rownames(y), taxa))
    z[, colnames(y)] <- y
    z
  })
  if (!is.null(taxonomy)) {
    taxonomy <- as.matrix(taxonomy)
    if (ncol(taxonomy) < 2L || is.null(rownames(taxonomy)) ||
        anyDuplicated(rownames(taxonomy)) || !all(taxa %in% rownames(taxonomy)))
      stop("taxonomy needs at least two ranks and unique row names covering all taxa.",
           call. = FALSE)
    if (!identical(colnames(taxonomy), paste0("Rank", seq_len(ncol(taxonomy)))))
      stop("taxonomy columns must be named Rank1, Rank2, ..., from higher to lower rank.",
           call. = FALSE)
    taxonomy <- taxonomy[taxa, , drop = FALSE]
    if (anyNA(taxonomy) || any(!nzchar(taxonomy)))
      stop("taxonomy cannot contain missing or empty rank labels; use explicit unclassified labels.",
           call. = FALSE)
  }
  list(counts = counts, X = X, X.index = X.index, cluster.id = cluster.id,
       taxonomy = taxonomy)
}
