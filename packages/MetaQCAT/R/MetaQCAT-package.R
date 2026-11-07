#' Wild-score bootstrap meta-analysis of microbiome associations
#'
#' MetaQCAT combines covariate-adjusted, null-based score summaries across
#' independent studies. Count-composition and observed-zero components support
#' fixed-effect and random-effect tests and Cauchy omnibus combinations.
#' Calibration uses a fixed-budget subject-level Rademacher multiplier bootstrap
#' with fixed covariance estimates. This is an asymptotic procedure, not exact
#' disease-label permutation. See [MetaQCAT()] and the package vignette for
#' input conventions, assumptions, diagnostics and validation scope.
#'
#' @name MetaQCAT-package
#' @importFrom Rcpp evalCpp
#' @useDynLib MetaQCAT, .registration = TRUE
"_PACKAGE"
