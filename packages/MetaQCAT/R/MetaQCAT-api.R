#' Cluster-score bootstrap for microbiome composition meta-analysis
#'
#' Combine study-specific composition and zero-pattern score tests using
#' fixed-effect, random-effect, or omnibus methods. Bootstrap weights are
#' sampled independently by subject cluster, retaining the observed relation
#' between the tested covariate and adjustment variables.
#'
#' @param counts A list of nonnegative integer-valued, sample-by-taxon count
#'   matrices, one per study. Taxon column names are required. Different study
#'   columns are aligned by name; absent taxa are treated as zero counts.
#' @param X A list of numeric covariate matrices. Rows must match `counts` in
#'   number and order. Do not include an intercept; one is added internally.
#'   Adjustment covariates may differ between studies. Named study lists are
#'   aligned by study name. Sample alignment is checked before all inputs are
#'   sorted together by stable study, sample and taxon names for reproducibility.
#' @param X.index Column indices in each `X` matrix for the covariate(s) to test,
#'   before adding an intercept. If omitted, all columns of the first matrix
#'   are tested. Specify this explicitly when adjustment covariates are present.
#' @param cluster.id A list of subject/cluster ID vectors, one per study, or
#'   `NULL` for independent samples. Entries may be `NULL` for individual
#'   independent studies. Each supplied vector must have one nonmissing ID
#'   per sample. IDs are converted to factors with unused levels removed.
#'   With `NULL`, unique sample row names serve as singleton subject IDs.
#'   Unnamed samples receive positional IDs; supply stable names to retain
#'   identical keyed draws when reordering the data.
#' @param taxonomy A taxon-by-rank matrix with unique row names covering all
#'   count-matrix taxa. Columns must be named `Rank1`, `Rank2`, etc., from
#'   higher to lower rank. Extra rows are ignored. With `NULL`, perform one
#'   composition-wide test; otherwise test parents with multiple child taxa.
#' @param method One of `"FE-MetaQCAT"`, `"RE-MetaQCAT"`, `"FE-MetaQCAT*"`,
#'   `"RE-MetaQCAT*"`, `"FE-MetaQCAT-O"`, or `"RE-MetaQCAT-O"`.
#'   Unstarred methods use composition information, starred methods use zero
#'   indicators, and `-O` combines component p-values using a Cauchy combination.
#' @param min.depth Remove samples with total count less than or equal to this
#'   threshold, within each study. Default: 0 (remove zero-depth samples).
#' @param min.prev Prevalence threshold between 0 and 1. For composition tests, keep
#'   taxa with nonzero proportion at least this value and some positive counts.
#'   For zero tests, require the zero proportion strictly between `min.prev`
#'   and `1 - min.prev`. Default: 0.1. Omnibus tests may use two values, for
#'   composition and zero components respectively. This is not a maximum
#'   zero-proportion threshold.
#' @param n.boot Positive integer number of cluster-level Rademacher score
#'   bootstrap replicates. All requested replicates are run; default: 999.
#' @param fdr.alpha Benjamini-Hochberg FDR threshold for lineage tests.
#'   Default: 0.05. Uses bootstrap p-values; failed lineage tests are not
#'   declared significant.
#' @param save.coef Whether to return descendant-taxon one-step score summaries
#'   and variances. Default: `FALSE`. The composition reference taxon is always
#'   reported when its test succeeds.
#' @param n.threads Maximum number of forked workers for lineage tests.
#'   Default: 1. Windows falls back to one worker with a warning.
#' @param seed Integer random seed for rarefaction and bootstrap weights.
#'   Default: 1234. One set of weights is shared across lineages and component
#'   tests, independently of worker scheduling.
#'
#' @details
#' The zero component is based on rarefied counts. After depth filtering, each
#' study is rarefied once to its own minimum sample depth, before aggregation
#' into lineages. This follows MetaPALM, not the older MetaQCAT implementation.
#' Rarefaction uses sequential hypergeometric sampling without expanding reads;
#' its distribution is unchanged, but exact seeded draws can differ from MetaPALM.
#' Cluster weights are generated after depth filtering and before component-
#' specific screening. A cluster's repeated samples share one weight per
#' replicate. Exposures and adjustment covariates are not permuted; users must
#' explicitly supply any within- or between-cluster exposure columns to test.
#' This is an asymptotic score bootstrap, not an exact randomization test.
#' The reduced multinomial and logistic null fits use unpenalized maximum
#' likelihood, with nuisance-score residuals recorded in the diagnostics.
#' The old bias-reduced composition fit is not used: accurate solution of the
#' ordinary nuisance equations is required for ordinary-score recovery.
#' Subject contributions are not centered, and the fixed covariance is their
#' unscaled cross-product. All requested draws must be finite.
#'
#' Omnibus results use the original equal-weight Cauchy combination of the
#' marginal bootstrap p-values, not a second joint rank calibration. An exact
#' component p-value of one gives a combined p-value of one. A zero gives zero;
#' the simultaneous zero/one asymptotic underflow case is handled by clipping
#' to `[1e-15, 1 - 1e-15]` before evaluating the tangent. Bootstrap marginal
#' p-values are positive because they use `(1 + exceedances)/(n.boot + 1)`.
#'
#' For the composition `"RE-MetaQCAT"` test, the reference-invariant
#' heterogeneity kernel is \eqn{B=(I_K+\mathbf{1}\mathbf{1}^T)\otimes I_q},
#' where \eqn{K} is the number of nonreference taxa and \eqn{q} is the
#' number of tested covariates.
#' One eligible taxon is selected as a common reference across studies and
#' retained for bootstrap calibration. Provided the same taxa are eligible and
#' required covariance matrices are positive definite, changing the reference leaves the RE
#' composition test unchanged up to numerical tolerance. In every method,
#' singular or non-positive-definite covariance or deficient pooled information
#' causes an explicit nonestimability failure rather than a pseudoinverse or
#' ridge. The zero-pattern RE component retains its identity kernel. Both RE
#' components use \eqn{\{\max(0,U_\sigma)\}^2/V_\sigma}, with the score and
#' information summed across studies before taking the positive part, in the
#' observed statistic and every bootstrap draw.
#'
#' Taxa absent from a study's count matrix are treated as zero, not as missing
#' measurements. Do not use this convention for unmeasured taxa. Studies must
#' be independent of one another; `cluster.id` only accounts for dependence
#' within each study. Overlap of subjects across studies is not accounted for.
#' Failure messages are retained rather than silently
#' removing lineage results. Finite p-values do not establish statistical
#' calibration; see the vignette for design and simulation considerations.
#'
#' @return Without `taxonomy`, component fields `pvalue.comp`, `pvalue.pa`,
#'   and/or `pvalue.comb`, named `Asy` and/or `Boot`; score statistics; and
#'   requested coefficient estimates. With `taxonomy`, a list containing:
#'   * `lineage.pval`: calibration-by-lineage p-value matrix;
#'   * `lineage.qval`: corresponding BH-adjusted p-values;
#'   * `lineage.coef`: detailed component results, including optional coefficients;
#'   * `sig.lineage`: BH-significant lineage names;
#'   * `failed.lineage`: named failure reasons. Failed tests retain `NA` columns.
#'   Successful RE composition results also include `reference.taxon` and
#'   `re.diagnostics` (kernel and coefficient ordering, score components,
#'   covariance rank/conditioning, seed, and requested/actual bootstrap counts).
#'   The `bootstrap` field records seed, requested/completed replicates, and
#'   independent cluster counts by study.
#'   The `score.summary.comp` and `score.summary.pa` fields retain subject-level
#'   `psi.lst` matrices, fixed covariances, coordinate maps, kernels, stable IDs,
#'   and numerical diagnostics. Successful tests have `status = "ok"`; failed
#'   tests have `status = "nonestimable"` and retain `NA` p-values.
#'   In lineage scans, inspect these fields within `lineage.coef[[name]]`.
#' @seealso [CRC_data], [IBD_data], [Score_test_QCAT_meta],
#'   `vignette("MetaQCAT", package = "MetaQCAT")`
#' @references
#' Tang ZZ, Chen G, Alekseyenko AV, Li H (2017). A general framework for
#' association analysis of microbial communities on a taxonomic tree.
#' Bioinformatics. \doi{10.1093/bioinformatics/btw804}.
#'
#' Lee S, Teslovich TM, Boehnke M, Lin X (2013). General framework for
#' meta-analysis of rare variants in sequencing association studies.
#' American Journal of Human Genetics. \doi{10.1016/j.ajhg.2013.05.010}.
#' @examples
#' data("CRC_data")
#' names(CRC_data$counts)
#' # See the vignette for runnable CRC and clustered IBD analyses.
#' @export
#' @import stats
MetaQCAT <- function(counts, X, X.index, cluster.id = NULL, taxonomy = NULL,
                     method = "FE-MetaQCAT", min.depth = 0, min.prev = 0.1,
                     n.boot = 999L, fdr.alpha = 0.05,
                     save.coef = FALSE, n.threads = 1, seed = 1234) {
  .check_arg_names(sys.call(), names(formals(MetaQCAT)))
  .check_scalar(min.depth, "min.depth")
  .check_scalar(fdr.alpha, "fdr.alpha", 0, 1)
  .check_scalar(n.threads, "n.threads", 1, .Machine$integer.max, integer = TRUE)
  .check_options(method, min.prev, n.boot, save.coef, seed)
  if (missing(X.index)) {
    if (!is.list(X) || !length(X)) stop("X must be a nonempty list.", call. = FALSE)
    X.index <- seq_len(ncol(as.matrix(X[[1L]])))
  }
  d <- .prepare_inputs(counts, X, X.index, cluster.id, taxonomy)
  for (i in seq_along(d$counts)) {
    keep <- rowSums(d$counts[[i]]) > min.depth
    d$counts[[i]] <- d$counts[[i]][keep, , drop = FALSE]
    d$X[[i]] <- d$X[[i]][keep, , drop = FALSE]
    if (!is.null(d$cluster.id) && !is.null(d$cluster.id[[i]]))
      d$cluster.id[[i]] <- factor(d$cluster.id[[i]][keep])
  }
  keep.study <- vapply(d$counts, nrow, integer(1)) > 0L
  if (!any(keep.study)) stop("No studies remain after depth filtering.", call. = FALSE)
  if (any(!keep.study)) warning("Excluded empty studies after depth filtering: ",
    paste(names(d$counts)[!keep.study], collapse = ", "), call. = FALSE)
  counts <- d$counts[keep.study]; X <- d$X[keep.study]
  cluster.id <- .complete_cluster_ids(counts,
    if (is.null(d$cluster.id)) NULL else d$cluster.id[keep.study])
  weights <- .make_cluster_weights(cluster.id, n.boot = n.boot, seed = seed)
  X <- lapply(X, function(x) cbind(`(Intercept)` = 1, x))
  X.index <- d$X.index + 1L
  zeros <- if (grepl("\\*$|-O$", method))
    .with_seed(as.integer((as.double(seed) + 1) %% .Machine$integer.max),
               lapply(counts, .rarefy_counts)) else NULL
  if (grepl("\\*$", method)) counts <- NULL
  if (is.null(taxonomy)) {
    if (!is.null(zeros)) zeros <- lapply(zeros, function(y) (y == 0) * 1)
    res <- .run_meta_test(counts, zeros, X, X.index, cluster.id, weights,
                          method, min.prev, n.boot, save.coef, seed)
    if (!is.null(res$error)) warning("Test failed; returning NA: ", res$error, call. = FALSE)
    return(res)
  }
  .scan_lineages(counts, zeros, X, X.index, cluster.id, weights, d$taxonomy,
                 method, min.prev, n.boot, fdr.alpha, save.coef, n.threads, seed)
}
