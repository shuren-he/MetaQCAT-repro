# MetaQCAT

**Fixed- and random-effects meta-analysis of microbiome composition associations.**

MetaQCAT combines covariate-adjusted study-level score summaries to test
compositional-count and observed-zero associations across independent studies.
It supports independent samples and repeated observations within subjects,
using a null-imposed, cluster-level wild score bootstrap.

## Installation

Install the latest package from GitHub:

```r
install.packages("remotes")  # Run once if needed
remotes::install_github("shuren-he/MetaQCAT")
library(MetaQCAT)
```

R >= 3.6.0 and a working C++ toolchain are required for source installation.
Rtools is needed on Windows; macOS users may need the Xcode Command Line Tools.

## Quick start

This small example tests a binary disease indicator while adjusting for age
in two independent studies. The low bootstrap budget is for demonstration.

```r
library(MetaQCAT)

set.seed(31)
counts <- X <- vector("list", 2L)
names(counts) <- names(X) <- c("Study_A", "Study_B")

for (s in seq_along(counts)) {
  counts[[s]] <- t(rmultinom(80L, size = 100L, prob = c(.4, .35, .25)))
  dimnames(counts[[s]]) <- list(
    paste0("sample_", seq_len(80L)),
    paste0("taxon_", seq_len(3L))
  )
  X[[s]] <- cbind(
    disease = rep(0:1, 40L),
    age = rnorm(80L, mean = 50, sd = 10)
  )
  rownames(X[[s]]) <- rownames(counts[[s]])
}

fit <- MetaQCAT(
  counts = counts,
  X = X,
  X.index = 1L,
  method = "FE-MetaQCAT",
  n.boot = 199L,
  seed = 42L
)

fit$pvalue.comp    # Asy and Boot p-values
fit$bootstrap      # Bootstrap settings and cluster information
```

## Available tests

| Method | Component | Calibration |
| --- | --- | --- |
| `FE-MetaQCAT` | Compositional-count | Asymptotic and bootstrap |
| `FE-MetaQCAT*` | Observed-zero | Asymptotic and bootstrap |
| `FE-MetaQCAT-O` | Omnibus combination of the two FE components | Asymptotic and bootstrap |
| `RE-MetaQCAT` | Compositional-count with a heterogeneity score | Bootstrap |
| `RE-MetaQCAT*` | Observed-zero with a heterogeneity score | Bootstrap |
| `RE-MetaQCAT-O` | Omnibus combination of the two RE components | Bootstrap |

The compositional-count component concerns child-taxon composition conditional
on the lineage total, not the lineage's total abundance in the microbiome.
The zero component uses zero indicators after within-study rarefaction; observed
zeros do not necessarily represent biological absence. The `-O` tests use an
equal-weight Cauchy combination of component p-values.

## Inputs and main arguments

- `counts`: one sample-by-taxon, nonnegative integer count matrix per study,
  supplied as a list. Use unique taxon column names and stable sample row names.
- `X`: one numeric design matrix per study, with rows matching `counts`.
  **Do not add an intercept.** Include the exposure and adjustment variables.
- `X.index`: column index or indices of the exposure(s) to test. Other columns
  are adjustment covariates. Specify this explicitly when covariates are present.
  Tested columns must have the same meaning and ordering across studies.
- `cluster.id`: one subject-ID vector per study for repeated samples.
  Repeated observations from the same subject share an ID. Use `NULL` for
  independent samples.
- `taxonomy`: optional taxon-by-rank matrix for lineage scans. Row names must
  cover the count-table taxa; columns must be `Rank1`, `Rank2`, etc., from
  higher to lower rank, with at least two ranks.
- `min.depth`: remove samples whose total count is **less than or equal to**
  this threshold. Default: `0`.
- `min.prev`: component-specific prevalence threshold. Default: `0.1`.
  Count tests retain taxa with positive-count prevalence at least this value;
  zero tests require zero prevalence strictly between `min.prev` and
  `1 - min.prev`. This is not a maximum zero-fraction parameter.
- `n.boot`: requested bootstrap draws, all of which are run. Default: `999L`.
- `fdr.alpha`: BH-adjusted significance threshold for lineage scans.
  Default: `0.05`.
- `n.threads`: workers for lineage scans. Default: `1`; Windows uses one worker.
- `seed`: random seed for bootstrap weights and rarefaction.
  Default: `1234`.

Taxa are aligned by name rather than column position. A taxon omitted from one
study's count table is treated as zero, **not** as an unmeasured feature.
Handle genuinely unmeasured taxa explicitly before using the package.

## Packaged CRC and IBD data

Two processed example datasets are included:

```r
data("CRC_data", package = "MetaQCAT")
data("IBD_data", package = "MetaQCAT")

names(CRC_data)
names(IBD_data)

?CRC_data
?IBD_data
```

Both contain `counts`, `metadata`, `subject_ids`, `taxonomy`, `disease`
and `adjust`. IBD also contains `design`, combining the disease indicator
with study-specific adjustment covariates.

### Family-level CRC scan

Subset and rename the taxonomy ranks to test genus composition within families:

```r
tax <- CRC_data$taxonomy[, c("Rank5", "Rank6"), drop = FALSE]
colnames(tax) <- c("Rank1", "Rank2")

crc_fit <- MetaQCAT(
  counts = CRC_data$counts,
  X = CRC_data$disease,
  X.index = 1L,
  taxonomy = tax,
  method = "FE-MetaQCAT-O",
  n.boot = 9999L,
  seed = 1234L
)

crc_fit$lineage.pval
crc_fit$lineage.qval
crc_fit$sig.lineage
crc_fit$failed.lineage
```

This CRC design tests disease status without additional adjustment covariates.

### Repeated-sample IBD analysis

Use the full adjusted design and the subject IDs:

```r
tax <- IBD_data$taxonomy[, c("Rank5", "Rank6"), drop = FALSE]
colnames(tax) <- c("Rank1", "Rank2")

ibd_fit <- MetaQCAT(
  counts = IBD_data$counts,
  X = IBD_data$design,
  X.index = 1L,
  cluster.id = IBD_data$subject_ids,
  taxonomy = tax,
  method = "RE-MetaQCAT-O",
  n.boot = 9999L,
  seed = 1234L
)

ibd_fit$lineage.qval["Bootstrap-RE-MetaQCAT-O", ]
```

These are usage examples, not commands that reproduce the paper's complete
analysis protocol.

## Bootstrap and interpretation

Labels and covariates are not permuted. Independent subject clusters receive
Rademacher multipliers; repeated samples within a cluster share a multiplier.
Weights are reused across components and lineages, while null fits and covariance
estimates remain fixed. Within- or between-subject exposures must be constructed
explicitly in `X`; there is no `permute.strata` option.

Marginal bootstrap p-values use `(1 + exceedances) / (n.boot + 1)`. The
minimum attainable value is `1 / (n.boot + 1)`; choose the budget with the
intended significance and multiple-testing thresholds in mind. Omnibus tests
combine the marginal p-values rather than performing a second joint bootstrap
calibration.

This is an asymptotic score bootstrap, not an exact label-permutation test.
Studies must be independent: `cluster.id` addresses dependence within studies,
not subjects shared across studies. Count RE reference invariance applies to
a fixed eligible taxon set under the required full-rank covariance conditions;
it does not imply invariance to changing the retained taxa.

Returned coefficient summaries are **one-step score summaries**, not full-model
effect estimates. Nonestimable lineage tests retain `NA` values and failure
reasons; inspect `failed.lineage` and the component diagnostics.

## Documentation and paper code

- Function help: `?MetaQCAT`.
- Detailed tutorial and methodological notes:
  [package vignette](vignettes/MetaQCAT.Rmd).
- To install and build the local HTML vignette:
  `remotes::install_github("shuren-he/MetaQCAT", build_vignettes = TRUE)`
  (requires knitr, rmarkdown and Pandoc), then run
  `vignette("MetaQCAT", package = "MetaQCAT")`.
- Simulation, real-data and figure workflows:
  [MetaQCAT-repro](https://github.com/shuren-he/MetaQCAT-repro).

## License

GPL (>= 3); see [DESCRIPTION](DESCRIPTION) and [LICENSE.md](LICENSE.md).
