# One frozen configuration source for opt-in validation, not R CMD check.
# Profiles only set budgets/designs; they do not launch work when sourced.
validation_config <- function(profile = "smoke") {
  stopifnot(profile %in% c("smoke", "principal", "pipeline", "power"))
  list(
    profile = profile, seed = 20261101L,
    package.version = as.character(utils::packageVersion("MetaQCAT")),
    dependency.versions = vapply(c("brglm2", "Rcpp", "RcppArmadillo", "parallelly"),
      function(p) as.character(utils::packageVersion(p)), character(1)),
    rng = c("Mersenne-Twister", "Inversion", "Rejection"),
    n.rep = if (profile == "smoke") 2L else 5000L,
    n.boot = if (profile == "smoke") 19L else 1999L,
    n.threads = 1L, alpha = c(.05, .01),
    subject.counts = c(40L, 80L, 120L),
    size.multipliers = c(application_like = 1L, larger = 4L),
    covariate.strata = 2L, disease.probability.independent = c(.5, .5),
    disease.probability.associated = c(.2, .8),
    repeat.sizes = 1:4, dirichlet.concentration = .5,
    fixed.depth = 40L, taxa = 4L,
    X.index = 1L, covariates = "saturated stratum indicator",
    pipeline = if (profile %in% c("pipeline", "power")) "full" else "fixed",
    min.depth = 0, min.prev = if (profile %in% c("pipeline", "power")) .05 else 0,
    rarefaction = profile %in% c("pipeline", "power"),
    rank.rule = "strict backend rank/conditioning checks; no pseudoinverse or ridge",
    scenarios = if (profile == "power") c("homogeneous", "one_flipped") else "null",
    signals = if (profile == "power") c("count", "zero") else "none",
    signal.log.multiplier = .6, signal.add.probability = .2,
    methods = c("FE-MetaQCAT", "RE-MetaQCAT", "FE-MetaQCAT*",
                "RE-MetaQCAT*", "FE-MetaQCAT-O", "RE-MetaQCAT-O"),
    multiplicity = "one lineage per simulated dataset; no FDR-control claim",
    matched.application.budget = 9999L,
    validation.status = "planned; execution smoke checks are not calibration evidence"
  )
}
