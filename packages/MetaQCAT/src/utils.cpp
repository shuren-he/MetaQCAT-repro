// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <sstream>

struct MatrixCheck { unsigned int rank; double condition; double tolerance; };

// Ordinary inverses only: unsupported rank/conditioning is nonestimability.
arma::mat solve_full(const arma::mat& a, const arma::mat& b,
                     const std::string& context, MatrixCheck& check) {
  if (a.n_rows == 0 || a.n_rows != a.n_cols || a.n_rows != b.n_rows ||
      !a.is_finite() || !b.is_finite())
    Rcpp::stop(context + " has invalid dimensions or non-finite entries.");
  arma::vec singular;
  if (!arma::svd(singular, a)) Rcpp::stop(context + " SVD failed.");
  check.tolerance = 64.0 * std::numeric_limits<double>::epsilon() * a.n_rows * singular.max();
  check.rank = arma::accu(singular > check.tolerance);
  check.condition = singular.min() > 0 ? singular.max() / singular.min() :
    std::numeric_limits<double>::infinity();
  if (check.rank != a.n_rows || !std::isfinite(check.condition)) {
    std::ostringstream msg;
    msg << context << " is not estimable (rank " << check.rank << "/" << a.n_rows
        << ", condition number " << check.condition << ", singular-value tolerance "
        << check.tolerance << ").";
    Rcpp::stop(msg.str());
  }
  arma::mat out;
  if (!arma::solve(out, a, b, arma::solve_opts::no_approx) || !out.is_finite())
    Rcpp::stop(context + " solve failed.");
  return out;
}

arma::mat solve_spd(const arma::mat& a, const arma::mat& b,
                    const std::string& context, MatrixCheck& check) {
  if (a.n_rows == 0 || a.n_rows != a.n_cols || a.n_rows != b.n_rows ||
      !a.is_finite() || !b.is_finite())
    Rcpp::stop(context + " has invalid dimensions or non-finite entries.");
  const double symmetry_tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
    a.n_rows * arma::norm(a, "inf");
  if (arma::norm(a - a.t(), "inf") > symmetry_tolerance)
    Rcpp::stop(context + " is not symmetric within numerical tolerance.");
  arma::mat symmetric = 0.5 * (a + a.t()); // roundoff only, after checking
  arma::vec eigenvalues;
  if (!arma::eig_sym(eigenvalues, symmetric)) Rcpp::stop(context + " eigendecomposition failed.");
  check.tolerance = 64.0 * std::numeric_limits<double>::epsilon() * a.n_rows *
    arma::max(arma::abs(eigenvalues));
  check.rank = arma::accu(eigenvalues > check.tolerance);
  check.condition = eigenvalues.min() > 0 ? eigenvalues.max() / eigenvalues.min() :
    std::numeric_limits<double>::infinity();
  if (check.rank != a.n_rows || !std::isfinite(check.condition)) {
    std::ostringstream msg;
    msg << context << " is not numerically positive definite / not estimable (rank "
        << check.rank << "/" << a.n_rows << ", condition number " << check.condition
        << ", eigenvalue tolerance " << check.tolerance << ").";
    Rcpp::stop(msg.str());
  }
  arma::mat out;
  if (!arma::solve(out, symmetric, b, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) ||
      !out.is_finite()) Rcpp::stop(context + " symmetric solve failed.");
  return out;
}

Rcpp::List matrix_diagnostics(const MatrixCheck& c) {
  return Rcpp::List::create(Rcpp::Named("rank") = c.rank,
    Rcpp::Named("condition.number") = c.condition, Rcpp::Named("tolerance") = c.tolerance);
}

arma::mat count_kernel(const arma::uword dimension, const unsigned int q) {
  if (q == 0 || dimension == 0 || dimension % q != 0)
    Rcpp::stop("Count coefficient dimension is not divisible by q.");
  const arma::uword k = dimension / q;
  return arma::kron(arma::eye(k, k) + arma::ones(k, k), arma::eye(q, q));
}

double re_qcat_statistic_value(const double u, const double v) {
  if (!std::isfinite(u) || !std::isfinite(v) || v <= 0)
    Rcpp::stop("Degenerate random-effect information.");
  return std::pow(std::max(0.0, u), 2) / v;
}

// [[Rcpp::export]]
Rcpp::List re_qcat_components(const arma::vec& beta_hat, const arma::mat& cov_beta,
                              const unsigned int q) {
  arma::mat kernel = count_kernel(beta_hat.n_elem, q);
  MatrixCheck check;
  arma::mat precision = solve_spd(cov_beta, arma::eye(beta_hat.n_elem, beta_hat.n_elem),
    "Random-effect sampling covariance", check);
  arma::vec z = precision * beta_hat;
  arma::mat wb = precision * kernel;
  const double u = 0.5 * (arma::as_scalar(z.t() * kernel * z) - arma::trace(wb));
  const double v = 0.5 * arma::trace(wb * wb);
  if (!std::isfinite(u) || !std::isfinite(v) || v <= 0)
    Rcpp::stop("Degenerate random-effect score information.");
  return Rcpp::List::create(Rcpp::Named("U_theta") = u, Rcpp::Named("V_theta") = v,
    Rcpp::Named("rank") = check.rank, Rcpp::Named("condition_number") = check.condition,
    Rcpp::Named("tolerance") = check.tolerance, Rcpp::Named("kernel") = kernel);
}

// [[Rcpp::export]]
double re_qcat_statistic(const double U_theta, const double V_theta) {
  return re_qcat_statistic_value(U_theta, V_theta);
}

// Common observed/bootstrap evaluator. Every matrix solve is outside the draw
// dimension. R supplies local-to-global maps after harmonizing ALR references.
Rcpp::List evaluate_meta(const Rcpp::List& beta_list, const Rcpp::List& cov_list,
                         const Rcpp::List& maps, const unsigned int dimension,
                         const unsigned int q, const std::string& method,
                         const std::string& component) {
  const unsigned int studies = beta_list.size();
  if (studies == 0 || cov_list.size() != studies || maps.size() != studies || dimension == 0 || q == 0)
    Rcpp::stop("Invalid study lists or coefficient dimensions.");
  const bool re = method.compare(0, 2, "RE") == 0;
  if (!re && method.compare(0, 2, "FE") != 0) Rcpp::stop("Unknown meta-analysis method.");
  const bool composition = component == "count";
  if (!composition && component != "zero" && component != "poisson") Rcpp::stop("Unknown score component.");
  const bool one_sided = composition || component == "zero";
  const arma::mat first = Rcpp::as<arma::mat>(beta_list[0]);
  const arma::uword draws = first.n_cols;
  if (draws == 0) Rcpp::stop("No observed or bootstrap summaries.");
  arma::mat information(dimension, dimension, arma::fill::zeros);
  arma::mat pooled_score(dimension, draws, arma::fill::zeros);
  arma::rowvec u(draws, arma::fill::zeros);
  double v = 0;
  Rcpp::List kernels(studies), diagnostics(studies);
  Rcpp::IntegerVector ranks(studies);
  Rcpp::NumericVector conditions(studies), tolerances(studies);
  for (unsigned int j = 0; j < studies; ++j) {
    arma::mat beta = Rcpp::as<arma::mat>(beta_list[j]);
    arma::mat covariance = Rcpp::as<arma::mat>(cov_list[j]);
    arma::uvec index = Rcpp::as<arma::uvec>(maps[j]);
    const arma::uword d = beta.n_rows;
    arma::uvec unique_index = arma::unique(index);
    if (d == 0 || beta.n_cols != draws || !beta.is_finite() || index.n_elem != d ||
        arma::any(index < 1) || arma::any(index > dimension) || unique_index.n_elem != index.n_elem)
      Rcpp::stop("Incompatible or non-finite study summaries / taxon-coordinate map.");
    index -= 1;
    MatrixCheck check;
    arma::mat precision = solve_spd(covariance, arma::eye(d, d),
      "Study " + std::to_string(j + 1) + " sampling covariance", check);
    diagnostics[j] = matrix_diagnostics(check);
    ranks[j] = check.rank; conditions[j] = check.condition; tolerances[j] = check.tolerance;
    information.submat(index, index) += precision;
    pooled_score.rows(index) += precision * beta;
    arma::mat kernel = composition ? count_kernel(d, q) : arma::eye(d, d);
    kernels[j] = kernel;
    if (re) {
      arma::mat wb = precision * kernel;
      arma::mat z = precision * beta;
      u += 0.5 * (arma::sum(z % (kernel * z), 0) - arma::trace(wb));
      v += 0.5 * arma::trace(wb * wb);
    }
  }
  MatrixCheck pooled_check;
  arma::mat pooled_covariance = solve_spd(information, arma::eye(dimension, dimension),
    "Pooled mean information", pooled_check);
  arma::mat estimates = pooled_covariance * pooled_score;
  arma::rowvec statistics = arma::sum(pooled_score % estimates, 0);
  if (re) {
    if (!std::isfinite(v) || v <= 0) Rcpp::stop("Degenerate random-effect score information.");
    // Count and zero components share the nonnegative-variance boundary.
    // Clip the summed score once, for the observed statistic and every draw.
    // The separate legacy Poisson helper retains its existing convention.
    arma::rowvec used_u = one_sided ? arma::clamp(u, 0.0, std::numeric_limits<double>::infinity()) : u;
    statistics += arma::square(used_u) / v;
  }
  if (!statistics.is_finite()) Rcpp::stop("Non-finite observed or bootstrap statistic.");
  return Rcpp::List::create(Rcpp::Named("statistics") = statistics,
    Rcpp::Named("pooled.estimates") = estimates, Rcpp::Named("pooled.covariance") = pooled_covariance,
    Rcpp::Named("re.score") = u, Rcpp::Named("re.information") = v,
    Rcpp::Named("kernel.lst") = kernels, Rcpp::Named("covariance.diagnostics") = diagnostics,
    Rcpp::Named("re.diagnostics") = Rcpp::List::create(
      Rcpp::Named("kernel") = composition ? "isotropic-contrast" : "identity",
      Rcpp::Named("coefficient.order") = "taxon-major", Rcpp::Named("q") = q,
      Rcpp::Named("U_theta") = u[0], Rcpp::Named("V_theta") = v,
      Rcpp::Named("positive.part.active") = one_sided && u[0] > 0,
      Rcpp::Named("variance.score.formula") = one_sided ? "positive-part" : "untruncated",
      Rcpp::Named("covariance.rank") = ranks, Rcpp::Named("covariance.condition.number") = conditions,
      Rcpp::Named("covariance.tolerance") = tolerances, Rcpp::Named("pooled.rank") = pooled_check.rank,
      Rcpp::Named("pooled.condition.number") = pooled_check.condition,
      Rcpp::Named("pooled.tolerance") = pooled_check.tolerance));
}

arma::uvec tested_indices(const arma::vec& target, const arma::uword p, const arma::uword taxa) {
  arma::vec unique_target = arma::unique(target);
  if (target.n_elem == 0 || !target.is_finite() || arma::any(target < 1) || arma::any(target > p) ||
      arma::any(target != arma::floor(target)) || unique_target.n_elem != target.n_elem)
    Rcpp::stop("Invalid tested coefficient indices.");
  arma::uvec out(taxa * target.n_elem);
  for (arma::uword k = 0; k < taxa; ++k)
    for (arma::uword h = 0; h < target.n_elem; ++h) out[k * target.n_elem + h] = k * p + target[h] - 1;
  return out;
}

struct SubjectSummary {
  arma::mat cluster_score, bread, sensitivity, psi;
  arma::vec beta;
  arma::mat covariance;
  Rcpp::List diagnostics;
};

SubjectSummary subject_summary(const arma::mat& score, const arma::mat& h,
                                const arma::uvec& tested, SEXP cluster, const std::string& context) {
  if (score.n_rows == 0 || score.n_cols != h.n_rows || h.n_rows != h.n_cols ||
      !score.is_finite() || !h.is_finite()) Rcpp::stop(context + " has invalid/non-finite score arrays.");
  arma::uvec mask(score.n_cols, arma::fill::ones); mask.elem(tested).zeros();
  arma::uvec nuisance = arma::find(mask == 1);
  arma::mat sensitivity = h.submat(tested, tested), adjusted = score.cols(tested);
  MatrixCheck nuisance_check = {0, 1.0, 0.0};
  double nuisance_max = 0, nuisance_relative = 0;
  if (nuisance.n_elem) {
    arma::mat hbg = h.submat(tested, nuisance), hgg = h.submat(nuisance, nuisance);
    // P=H_beta,gamma H_gamma,gamma^{-1}, with correct nonsymmetric orientation.
    arma::mat projection = solve_full(hgg.t(), hbg.t(), context + " nuisance sensitivity", nuisance_check).t();
    sensitivity -= projection * h.submat(nuisance, tested);
    adjusted -= score.cols(nuisance) * projection.t();
    nuisance_max = arma::max(arma::abs(arma::sum(score.cols(nuisance), 0)));
    nuisance_relative = nuisance_max / std::max(1.0, arma::accu(arma::abs(score.cols(nuisance))));
  }
  MatrixCheck sensitivity_check;
  arma::mat bread = solve_full(sensitivity, arma::eye(tested.n_elem, tested.n_elem),
    context + " adjusted sensitivity", sensitivity_check);
  arma::mat cluster_scores;
  if (Rf_isNull(cluster)) cluster_scores = adjusted.t();
  else {
    Rcpp::IntegerVector ids(cluster);
    if (ids.size() != static_cast<int>(score.n_rows)) Rcpp::stop(context + " cluster IDs have wrong length.");
    int groups = 0;
    for (int id : ids) {
      if (id == NA_INTEGER || id < 1) Rcpp::stop(context + " cluster IDs must be nonmissing positive factor codes.");
      groups = std::max(groups, id);
    }
    arma::uvec seen(groups, arma::fill::zeros);
    cluster_scores.zeros(tested.n_elem, groups);
    for (arma::uword i = 0; i < score.n_rows; ++i) {
      seen[ids[i] - 1] = 1; cluster_scores.col(ids[i] - 1) += adjusted.row(i).t();
    }
    if (arma::any(seen == 0)) Rcpp::stop(context + " cluster IDs contain unused levels.");
  }
  // No centering or finite-sample scaling. Rows of Psi represent subjects.
  arma::mat psi = (bread * cluster_scores).t();
  arma::vec beta = arma::sum(psi, 0).t();
  arma::mat covariance = psi.t() * psi; // bread * meat * bread.t()
  if (!psi.is_finite() || !covariance.is_finite()) Rcpp::stop(context + " produced non-finite subject contributions.");
  Rcpp::List diagnostics = Rcpp::List::create(Rcpp::Named("n.subjects") = psi.n_rows,
    Rcpp::Named("nuisance.score.max") = nuisance_max, Rcpp::Named("nuisance.score.relative") = nuisance_relative,
    Rcpp::Named("nuisance.sensitivity") = matrix_diagnostics(nuisance_check),
    Rcpp::Named("adjusted.sensitivity") = matrix_diagnostics(sensitivity_check));
  return {cluster_scores, bread, sensitivity, psi, beta, covariance, diagnostics};
}

Rcpp::List pack_summary(const Rcpp::List& beta, const Rcpp::List& covariance,
                        const Rcpp::List& clusters, const Rcpp::List& bread, const Rcpp::List& psi,
                        const Rcpp::List& sensitivity, const Rcpp::List& diagnostics,
                        const Rcpp::List& maps, const unsigned int dimension, const unsigned int q,
                        const std::string& method, const std::string& component) {
  Rcpp::List meta = evaluate_meta(beta, covariance, maps, dimension, q, method, component);
  arma::rowvec statistics = Rcpp::as<arma::rowvec>(meta["statistics"]);
  arma::mat pooled = Rcpp::as<arma::mat>(meta["pooled.estimates"]), cov = Rcpp::as<arma::mat>(meta["pooled.covariance"]);
  arma::vec pooled_beta = pooled.col(0), pooled_variance = cov.diag();
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("score.stat.meta") = statistics[0],
    Rcpp::Named("cov.lst") = covariance, Rcpp::Named("cluster.score.lst") = clusters,
    Rcpp::Named("bread.lst") = bread, Rcpp::Named("psi.lst") = psi,
    Rcpp::Named("sensitivity.lst") = sensitivity, Rcpp::Named("score.diagnostics") = diagnostics,
    Rcpp::Named("covariance.diagnostics") = meta["covariance.diagnostics"],
    Rcpp::Named("kernel.lst") = meta["kernel.lst"], Rcpp::Named("nominal.df") = dimension,
    Rcpp::Named("asymptotic.df") = dimension);
  if (component == "zero") {
    out["alpha.hat"] = pooled_beta; out["var.alpha.hat"] = pooled_variance; out["est.alpha.lst"] = beta;
  } else {
    out["beta.hat"] = pooled_beta; out["var.beta.hat"] = pooled_variance; out["est.beta.lst"] = beta;
  }
  if (method.compare(0, 2, "RE") == 0) out["re.diagnostics"] = meta["re.diagnostics"];
  else out["score.pvalue"] = R::pchisq(statistics[0], dimension, 0, 0);
  return out;
}

// [[Rcpp::export]]
Rcpp::List score_bootstrap_meta(const Rcpp::List& component_summary, const Rcpp::List& weight_list,
                                const Rcpp::List& col_index_list, const unsigned int n_par_interest_meta,
                                const unsigned int q, const std::string& meta_method, const std::string& component) {
  Rcpp::List scores = component_summary["cluster.score.lst"], bread = component_summary["bread.lst"],
    covariance = component_summary["cov.lst"];
  const unsigned int studies = scores.size();
  if (studies == 0 || weight_list.size() != studies || bread.size() != studies)
    Rcpp::stop("Bootstrap score summary has inconsistent study counts.");
  arma::mat first_weights = Rcpp::as<arma::mat>(weight_list[0]);
  const arma::uword n_boot = first_weights.n_cols;
  if (n_boot == 0) Rcpp::stop("At least one bootstrap draw is required.");
  Rcpp::List beta_reps(studies);
  for (unsigned int j = 0; j < studies; ++j) {
    arma::mat s = Rcpp::as<arma::mat>(scores[j]), jinv = Rcpp::as<arma::mat>(bread[j]),
      weights = Rcpp::as<arma::mat>(weight_list[j]);
    if (s.n_rows == 0 || s.n_cols == 0 || jinv.n_rows != s.n_rows || jinv.n_cols != s.n_rows ||
        weights.n_rows != s.n_cols || weights.n_cols != n_boot || !s.is_finite() ||
        !jinv.is_finite() || !weights.is_finite()) Rcpp::stop("Incompatible/non-finite bootstrap score dimensions.");
    if (arma::any(arma::vectorise(arma::abs(weights) != 1.0))) Rcpp::stop("Multiplier weights must be Rademacher values (-1 or +1).");
    arma::mat psi_t = jinv * s, draws(s.n_rows, n_boot + 1);
    draws.col(0) = arma::sum(psi_t, 1); draws.cols(1, n_boot) = psi_t * weights;
    beta_reps[j] = draws;
  }
  Rcpp::List meta = evaluate_meta(beta_reps, covariance, col_index_list, n_par_interest_meta, q, meta_method, component);
  arma::rowvec statistics = Rcpp::as<arma::rowvec>(meta["statistics"]);
  Rcpp::NumericVector bootstrap_statistics(n_boot);
  for (arma::uword b = 0; b < n_boot; ++b) bootstrap_statistics[b] = statistics[b + 1];
  return Rcpp::List::create(Rcpp::Named("observed.statistic") = statistics[0],
    Rcpp::Named("bootstrap.statistic") = bootstrap_statistics, Rcpp::Named("study.beta.reps") = beta_reps,
    Rcpp::Named("re.score") = meta["re.score"], Rcpp::Named("re.information") = meta["re.information"]);
}

void count_score(const arma::vec& coefficients, const arma::mat& y, const arma::mat& x,
                  arma::mat& scores, arma::mat& h) {
  const arma::uword n = y.n_rows, k = y.n_cols - 1, p = x.n_cols;
  for (arma::uword i = 0; i < n; ++i) {
    arma::vec eta(k + 1, arma::fill::zeros);
    for (arma::uword a = 0; a < k; ++a) eta[a] = arma::dot(coefficients.subvec(a * p, (a + 1) * p - 1), x.row(i));
    eta -= eta.max();
    arma::vec probability = arma::exp(eta); probability /= arma::sum(probability);
    arma::vec means = probability.head(k);
    const double depth = arma::accu(y.row(i));
    scores.row(i) = arma::kron(y.row(i).head(k) - depth * means.t(), x.row(i));
    h += depth * arma::kron(arma::diagmat(means) - means * means.t(), x.row(i).t() * x.row(i));
  }
}

void zero_score(const arma::vec& coefficients, const arma::mat& y, const arma::mat& x,
                 arma::mat& scores, arma::mat& h) {
  const arma::uword n = y.n_rows, k = y.n_cols, p = x.n_cols;
  for (arma::uword i = 0; i < n; ++i) {
    arma::vec probability(k);
    for (arma::uword a = 0; a < k; ++a) {
      const double eta = arma::dot(coefficients.subvec(a * p, (a + 1) * p - 1), x.row(i));
      probability[a] = eta >= 0 ? 1.0 / (1.0 + std::exp(-eta)) : std::exp(eta) / (1.0 + std::exp(eta));
    }
    scores.row(i) = arma::kron(y.row(i) - probability.t(), x.row(i));
    h += arma::kron(arma::diagmat(probability % (1.0 - probability)), x.row(i).t() * x.row(i));
  }
}

Rcpp::List glm_summary(const Rcpp::List& x_list, const Rcpp::Nullable<Rcpp::List>& cluster_list,
                       const Rcpp::List& maps, const Rcpp::List& y_list, const Rcpp::List& coefficients,
                       const arma::vec& target, const unsigned int dimension, const std::string& method,
                       const bool zero) {
  const unsigned int studies = x_list.size();
  if (studies == 0 || y_list.size() != studies || coefficients.size() != studies) Rcpp::stop("Incompatible null-model study lists.");
  Rcpp::List ids;
  if (cluster_list.isNotNull()) { ids = cluster_list.get(); if (ids.size() != studies) Rcpp::stop("Wrong number of cluster-ID lists."); }
  Rcpp::List beta(studies), covariance(studies), clusters(studies), bread(studies), psi(studies), sensitivity(studies), diagnostics(studies);
  for (unsigned int j = 0; j < studies; ++j) {
    arma::mat x = Rcpp::as<arma::mat>(x_list[j]), y = Rcpp::as<arma::mat>(y_list[j]);
    if (x.n_rows != y.n_rows || y.n_cols < (zero ? 1u : 2u) || !x.is_finite() || !y.is_finite()) Rcpp::stop("Invalid count/design matrices.");
    const arma::uword k = y.n_cols - (zero ? 0 : 1), d = k * x.n_cols;
    arma::uvec tested = tested_indices(target, x.n_cols, k), mask(d, arma::fill::ones); mask.elem(tested).zeros();
    arma::uvec nuisance = arma::find(mask == 1);
    arma::vec coef = Rcpp::as<arma::vec>(coefficients[j]);
    if (coef.n_elem != nuisance.n_elem || !coef.is_finite()) Rcpp::stop("Non-finite or incorrectly ordered null coefficients.");
    arma::vec full_coef(d, arma::fill::zeros); full_coef.elem(nuisance) = coef;
    arma::mat scores(x.n_rows, d, arma::fill::zeros), h(d, d, arma::fill::zeros);
    if (zero) zero_score(full_coef, y, x, scores, h); else count_score(full_coef, y, x, scores, h);
    SEXP cluster = cluster_list.isNotNull() ? static_cast<SEXP>(ids[j]) : R_NilValue;
    SubjectSummary s = subject_summary(scores, h, tested, cluster, "Study " + std::to_string(j + 1));
    beta[j] = s.beta; covariance[j] = s.covariance; clusters[j] = s.cluster_score; bread[j] = s.bread;
    psi[j] = s.psi; sensitivity[j] = s.sensitivity; diagnostics[j] = s.diagnostics;
  }
  return pack_summary(beta, covariance, clusters, bread, psi, sensitivity, diagnostics, maps,
    dimension, target.n_elem, method, zero ? "zero" : "count");
}

// [[Rcpp::export]]
Rcpp::List score_test_stat_QCAT_meta(const Rcpp::List& X_list, const Rcpp::Nullable<Rcpp::List>& cluster_list,
                                     const Rcpp::List& col_index_list, const Rcpp::List& Y_list,
                                     const Rcpp::List& coef_List, const arma::vec& X_par_index,
                                     const unsigned int n_par_interest_beta, const std::string& meta_method) {
  return glm_summary(X_list, cluster_list, col_index_list, Y_list, coef_List, X_par_index, n_par_interest_beta, meta_method, false);
}

// [[Rcpp::export]]
Rcpp::List score_test_stat_zero_meta(const Rcpp::List& Z_list, const Rcpp::Nullable<Rcpp::List>& cluster_list,
                                     const Rcpp::List& col_index_list, const Rcpp::List& Y_bin_list,
                                     const Rcpp::List& coef_List, const arma::vec& Z_par_index,
                                     const unsigned int n_par_interest_alpha, const std::string& meta_method) {
  return glm_summary(Z_list, cluster_list, col_index_list, Y_bin_list, coef_List, Z_par_index, n_par_interest_alpha, meta_method, true);
}

// Poisson marginal-score variant retained for the existing lower-level API.
// [[Rcpp::export]]
Rcpp::List score_test_stat_meta(const Rcpp::List& X_list, const Rcpp::Nullable<Rcpp::List>& cluster_list,
                                const Rcpp::List& col_index_list, const Rcpp::List& Y_R_list,
                                const Rcpp::List& Y_I_list, const arma::uvec& X_par_index,
                                const unsigned int n_par_interest_beta, const std::string& meta_method) {
  const unsigned int studies = X_list.size();
  if (studies == 0 || Y_R_list.size() != studies || Y_I_list.size() != studies) Rcpp::stop("Incompatible Poisson-score study lists.");
  Rcpp::List ids;
  if (cluster_list.isNotNull()) { ids = cluster_list.get(); if (ids.size() != studies) Rcpp::stop("Wrong number of cluster-ID lists."); }
  Rcpp::List beta(studies), covariance(studies), clusters(studies), bread(studies), psi(studies), sensitivity(studies), diagnostics(studies);
  arma::vec target = arma::conv_to<arma::vec>::from(X_par_index);
  for (unsigned int j = 0; j < studies; ++j) {
    arma::mat x = Rcpp::as<arma::mat>(X_list[j]), residual = Rcpp::as<arma::mat>(Y_R_list[j]), information = Rcpp::as<arma::mat>(Y_I_list[j]);
    const arma::uword n = x.n_rows, k = residual.n_cols, p = x.n_cols;
    if (n == 0 || k == 0 || residual.n_rows != n || information.n_rows != n || information.n_cols != k ||
        !x.is_finite() || !residual.is_finite() || !information.is_finite()) Rcpp::stop("Invalid Poisson-score arrays.");
    arma::mat scores(n, k * p), h(k * p, k * p, arma::fill::zeros);
    for (arma::uword i = 0; i < n; ++i) scores.row(i) = arma::kron(residual.row(i), x.row(i));
    for (arma::uword a = 0; a < k; ++a) {
      arma::mat weighted = x.each_col() % information.col(a);
      h.submat(a * p, a * p, (a + 1) * p - 1, (a + 1) * p - 1) = x.t() * weighted;
    }
    SEXP cluster = cluster_list.isNotNull() ? static_cast<SEXP>(ids[j]) : R_NilValue;
    SubjectSummary s = subject_summary(scores, h, tested_indices(target, p, k), cluster, "Study " + std::to_string(j + 1));
    beta[j] = s.beta; covariance[j] = s.covariance; clusters[j] = s.cluster_score; bread[j] = s.bread;
    psi[j] = s.psi; sensitivity[j] = s.sensitivity; diagnostics[j] = s.diagnostics;
  }
  return pack_summary(beta, covariance, clusters, bread, psi, sensitivity, diagnostics, col_index_list,
    n_par_interest_beta, target.n_elem, meta_method, "poisson");
}
