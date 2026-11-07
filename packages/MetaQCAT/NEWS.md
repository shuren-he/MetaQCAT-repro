# MetaQCAT 0.2.1

- Unified the count and observed-zero RE boundary conventions: both use the
  positive part of the summed variance-component score in the observed
  statistic and every bootstrap draw. Their respective kernels are unchanged.
- RE-MetaQCAT* and RE-MetaQCAT-O p-values can differ from version 0.2.0.
  Previously saved results for these methods, including simulation summaries,
  must be recalculated using the updated statistic. FE methods and the
  compositional RE statistic are unchanged.

# MetaQCAT 0.2.0

- Replaced public label permutation with a null-imposed, subject-level
  Rademacher wild-score bootstrap. The interface uses `n.boot`; it has no
  `n.perm` or within/between permutation option. All draws are completed.
- Retained nuisance-adjusted subject contributions without centering. Observed
  and bootstrap summaries use the same contributions and fixed unscaled
  cluster-sandwich covariance. Strict rank failures are reported explicitly.
- Shared study/subject-keyed weights across components, FE/RE methods, lineages
  and worker schedules; independent observations are singleton subjects.
- Retained the induced ALR count kernel and the positive part of the aggregate
  count variance score. The zero component retains its identity kernel and
  untruncated legacy variance-score formula.
- Combined marginal plus-one bootstrap p-values using the existing equal-weight
  Cauchy rule, rather than the earlier experimental joint-rank omnibus.
- Added reproducibility metadata, deterministic tests, updated documentation,
  and opt-in calibration/power configurations. Large simulations and final
  application reruns are separate validation tasks, not completed by package
  checks. No broad calibration or FDR-control claim follows from execution tests.
- Preserved the pre-bootstrap source externally in `packages/MetaQCAT-old`.
  The archive keeps the original valid R package name internally.

# MetaQCAT 0.0.4.0000

- Created a standalone MetaQCAT package from miMeta 0.0.3.0000, with its own
  package metadata, native registration, documentation, vignette, and tests.
- Made the RE composition test invariant to the eligible reference taxon using
  the contrast-covariant kernel `B = (I_K + 1_K 1_K^T) ⊗ I_q` and the positive part
  of the total heterogeneity score. A common eligible reference is used across
  studies and resamples; the RE zero-pattern component is unchanged.
- Added `re.diagnostics` for the composition RE kernel, score, covariance
  conditioning, and resampling count. Singular/non-positive-definite covariance
  now fails explicitly instead of using a silent pseudoinverse. This change
  does not resolve conditional exchangeability of covariate-adjusted label
  permutations.

# miMeta 0.0.3.0000 (source package)

- Renamed the package from MetaQCAT to miMeta. The exported `MetaQCAT()`
  function and FE/RE-MetaQCAT method names remain unchanged for compatibility.

# MetaQCAT 0.0.3.0000 (prior package name)

## MetaPALM-based release

- Renamed the package, shared-library registration, generated Rcpp bindings,
  project file, and references to MetaQCAT. MetaPALM's original directory is
  retained outside this release as the reference implementation.
- Unified exported argument names: counts, X, X.index, cluster.id, taxonomy,
  method, min.prev, save.coef, and n.threads. Zero-indicator APIs use zeros.
  Old aliases and partial argument names are not accepted. Six low-level
  Score_test functions remain available with explicit intercept conventions.
- Added a public seed argument for rarefaction and resampling.

## Correctness and diagnostics

- Internal component dispatch uses named arguments, preventing accidental
  binding of taxon names or cluster IDs to prevalence parameters.
- Restored random-effect method labels in the retained Poisson backend.
- Aligned taxa and named study lists before analysis; checked sample alignment,
  covariates, taxonomy, tested indices, and cluster IDs before calling C++.
- Corrected empty-study detection after depth filtering and no-available-study
  handling in the composition component.
- Kept failed lineage columns as NA with explicit failure reasons rather than
  silently dropping them. Added lineage.pval and lineage.qval matrices and
  retained component results and reference taxa.
- Required an explicit clustered permutation scheme for FE as well as RE;
  checked exposure constancy for between-cluster and nontrivial within-cluster
  permutations. Removed unused factor levels before entering C++.
- Capped adaptive batches at the requested budget so remaining draws are not
  omitted. Rejected non-finite permutation statistics and degenerate RE
  covariance rather than returning a misleading numeric p-value.
- Replaced read-expansion rarefaction with sequential hypergeometric sampling:
  the sampling distribution is unchanged, but exact seeded draws may differ.
- Removed machine-specific compiler include paths and a deprecated finite check.

## Documentation and validation

- Added a styled HTML vignette with runnable CRC and clustered IBD examples,
  all six methods, result interpretation, failure diagnostics, and migration.
- Integrated current CRC/IBD dataset documentation into R/data_doc.R.
- Added automated regression tests for interfaces, all six independent and
  correlated resampling paths, dimensions, alignment, failures, and reproducibility.

## Scope of validation

The zero resampling interface in this MetaPALM-based source is connected; the
missing-interface issue in the older MetaQCAT source does not apply here.
Execution tests are not proofs of type-I error control. Cross-study subject
overlap, covariate-dependent exchangeability, sparse-data calibration, and
adaptive-resampling precision still require scientific sensitivity assessment.
