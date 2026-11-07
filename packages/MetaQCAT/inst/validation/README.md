# Opt-in wild-score bootstrap validation

These programs are not run during `R CMD check`. Install **MetaQCAT 0.2.1**
first, and deliberately choose a profile. Nothing in a successful software
check establishes type I error or FDR control.

```sh
Rscript run_simulation.R config.R smoke /absolute/path/to/smoke
Rscript run_simulation.R config.R principal /absolute/path/to/principal
Rscript run_simulation.R config.R pipeline /absolute/path/to/pipeline
Rscript run_simulation.R config.R power /absolute/path/to/power
Rscript benchmark_ibd.R --family=f__Enterobacteriaceae --boot=1000 \
  --seed=20261101 --output=/absolute/path/to/ibd-benchmark
```

`config.R` freezes seeds, bootstrap and dataset budgets, subject counts,
covariates, depth, filtering, rank policy, runtime settings and multiplicity
scope. Each result RDS additionally stores the actual installed package version,
path, dependency environment and dataset/bootstrap seeds. Archive the exact
package source and script revision with these objects before a publication run.
The scripts do not upload data, archive a release or update manuscript results.

The smoke profile requests two datasets and 19 draws per setting and checks
execution only. **Principal, pipeline and power profiles each request 5,000
datasets per setting and 1,999 draws**, and can be expensive. They are not
launched automatically. Use a new frozen configuration for budget changes;
do not increase a budget selectively to obtain significance.

## Correctly specified null benchmark

The principal profile crosses independent observations or unequal-size repeated
clusters with disease independent of or associated with a measured binary
stratum. Conditional disease probabilities are 0.5/0.5 or 0.2/0.8. Two subject
count scales assess convergence; the smaller first study has 40 subjects, so
one stratum/imbalance setting can have roughly eight controls. Actual generated
sample and independent-subject counts are retained.

Each subject has a latent Dirichlet composition centered on a study/stratum
mean. Counts are multinomial at fixed depth; repeated visits share the latent
composition. Unequal visit counts are drawn independently of disease and
latent composition within stratum. Disease has no conditional outcome effect.
The nuisance model is saturated in stratum, so both count means and taxon-zero
probabilities are representable. The benchmark fixes four eligible taxa and
the count reference, uses no extra rarefaction, and reports an unsupported
invariant feature as a failure rather than silently changing this design.

Count and zero null summaries are fitted once per dataset and reused for FE/RE
evaluation. All six tests use the same generated dataset and keyed multipliers.
Omnibus p-values are direct equal-weight Cauchy combinations of marginal
plus-one bootstrap p-values, not joint empirical-rank bootstrap combinations.

## Full pipeline and alternatives

The pipeline profile is a stylized sparse-count sensitivity design with
unequal study sizes/depths and the complete public filtering/rarefaction path.
It is not claimed to reproduce every CRC/IBD covariate or abundance pattern.
Preprocessing is performed once per public method fit and remains fixed over
its bootstrap draws. All methods share the same preprocessing seed and data.
Application-specific matching should be frozen separately before final reruns.

The power profile adds homogeneous or one-study-flipped signals with a random
baseline direction and target taxon, using paired data for every method.
Count-driven alternatives multiply the latent proportion; zero-driven
alternatives add a small binomial count before restoring total depth. The
latter also changes composition and is not called a pure zero-only alternative.
Interpret power only after null calibration is satisfactory.

## Outputs and interpretation

`replicates.csv` retains all six bootstrap p-values per dataset, failures,
requested/completed draws, subject/sample counts, seeds and paired-dataset
runtime. `summary.csv` reports rejection rates at 0.05 and 0.01, exact binomial
Monte Carlo intervals, MCSE, attempted/successful/failed replicates, and failure
rates. Conditional-on-success and all-attempted denominators are distinguished;
failed analyses are not silently discarded. Runtime is the complete paired
six-method dataset cost, not an independently measured method-specific cost.
`validation.rds` retains these tables, configuration and session information;
`failures.csv` preserves reasons. Inspect full p-value distributions, not only
threshold summaries. A confidence interval containing the nominal level is
not proof of calibration; systematic deviations need investigation at larger
sample sizes. Calibration of the old permutation backend is not a reference
truth. A paired old/new comparison requires frozen identical observed
statistics, data and preprocessing in the separate archived implementation.

The benchmark script times null fitting and fixed-fit bootstrap evaluation
separately on a real packaged IBD lineage, retaining explicit failures and cost
per 1,000 draws. It does not replace a full CRC/IBD rerun. Use matched fixed
budgets (the configuration records 9,999 as an initial application budget) for
full data and cohort omissions; assess the p-value floor against BH thresholds.
Global- and mixed-null full-lineage omnibus-plus-BH experiments are still
required before claiming FDR control. No such evidence is implied by this
single-lineage validation runner.
