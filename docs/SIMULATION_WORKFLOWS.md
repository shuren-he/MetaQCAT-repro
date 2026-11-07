# Paper simulation workflows

Run from the repository root after installing the bundled MetaQCAT 0.2.1 and
comparator dependencies (HMP, miLineage and MMUPHin). There is one generation
and preprocessing workflow, not separate generators for our tests and baselines.

## External inputs

Restore these files under `data/simu/`:

- mod_gdm_lst_Enterobacteriaceae.RData
- mod_gdm_lst_Peptostreptococcaceae.RData
- meta.rda, containing meta
- count.genus.rda, containing the aligned count.genus table

The model files contain five study-specific fitted generalized Dirichlet
parameter sets and taxon-order mappings. These inputs are not included in this
code-only release. No original GD-fitting script was located; generation is
reproducible conditional on the stored fitted parameters.

Study order is AUS, CHI, FRA, GER, USA; planned sample sizes are
109, 127, 114, 120, 104, respectively. These are generation-stage sample sizes,
not post-filtering real-data sample sizes. Enterobacteriaceae uses four genera:
Escherichia/Shigella, Citrobacter, Klebsiella and Enterobacter (excluding Hafnia).
Peptostreptococcaceae uses Peptostreptococcus, Intestinibacter and the unknown
Peptostreptococcaceae genus. Exact source column names are specified in the
entrypoints and validated against the count table.

## Run the paper analyses

```sh
Rscript scripts/simulation/CRC_Entero_Sim.R
Rscript scripts/simulation/CRC_Pepto_Sim.R
Rscript scripts/simulation/CRC_Entero_Type1.R
Rscript scripts/simulation/CRC_Pepto_Type1.R
```

Each replicate generates one dataset, then evaluates all 18 methods:

- FE/RE-MetaQCAT, FE/RE-MetaQCAT* and FE/RE-MetaQCAT-O;
- P, CD and RD, each with QCAT1, QCAT2, DM and ZIGDM.

Every method receives the same generated counts, labels, retained sample rows
and taxon identities. P combines single-study p-values using Fisher's method.
CD pools data after MMUPHin batch adjustment, preserving the simulated
case/control variable. RD pools the generated data directly. Native
method-specific rarefaction and internal eligibility rules remain part of
each test; QCAT2 uses rarefied data rather than replacing raw-count inputs for
all methods.

The four entrypoints set family constants and target, then source
simulation_workflow.R, whose settings, replicate loop, saving and summary code
are directly executable. simulation_helpers.R contains the generator and the
statistical calculations. Only the current MetaQCAT wild-bootstrap backend is
used.

## Default settings

| Setting | Default |
| --- | --- |
| Data / bootstrap seed | 2333 + replicate index |
| Power replicates per setting | 1,000 |
| Type I error replicates per design | 5,000 |
| MetaQCAT bootstrap draws | 300; no early stopping |
| Comparator permutation budget per fit | 300; RD-DM uses 1,000 |
| Requested workers | 7; capped by available cores, 1 on Windows |
| Significance level | 0.05 |
| MetaQCAT minimum prevalence | 0.1 |

Power has two mechanisms (proportional and count-addition) crossed with two
patterns (homogeneous and one-study-reversed). One taxon is randomly chosen.
Proportions are multiplied by 1 + Delta, where Delta is 0.75 for
Peptostreptococcaceae and 1.5 for Enterobacteriaceae. The initial direction is
chosen with probability 1/2; the heterogeneous setting reverses one study.
Count addition draws a fresh Binomial(trials, probability) value per affected
sample, with probabilities 0.05 and 0.005, respectively, and trials equal to
the rounded overall observed family mean depth. Modified count rows are
renormalized to their original simulated depths.

The null baseline draws proportions from the intercept-only study-specific
GD model. Each sample depth is Poisson with mean equal to that study's observed
family mean depth: neither a fixed depth nor a resampled observed sample depth.
Samples below 1.5 times the number of taxa are removed before any test.
Multinomial sampling generates counts from the proportions and depth.
Independent labels split retained samples as evenly as possible, assigning an
extra sample to cases when necessary.

Both independent and correlated type I error designs run by default, with
signal size zero. Correlated samples are randomly paired within study. Both
members' independent GD vectors are averaged with a third, shared GD baseline
vector. Labels are assigned to clusters. An odd remaining sample stays in a
singleton cluster. MetaQCAT uses these subject IDs for the wild bootstrap.
The comparator implementations do not accept cluster IDs; the correlated
experiment assesses their behavior without silently replacing them by
different cluster-aware tests.

For MetaQCAT, each count/zero null model is fitted once per study and replicate.
FE/RE reuse the scores and the count/zero components share keyed cluster
weights. Marginal p-values are (1 + exceedances)/(n.boot + 1); -O uses the
current package's Cauchy combination of these marginal p-values. RE count and
zero components both use the positive part of the summed variance-component
score. DM preserves its exceedances/successful-permutations rule without a
plus-one correction, and miLineage retains its own permutation backend.
The configured budgets are therefore not a single common resampling scheme.

## Small checks and controls

Environment controls are SIM_REPS, SIM_BOOTS, SIM_PERMS, SIM_RD_DM_PERMS,
SIM_THREADS and SIM_OUTPUT. They apply to power and type I error entrypoints.

```sh
SIM_REPS=2 SIM_BOOTS=19 SIM_PERMS=19 SIM_RD_DM_PERMS=19 SIM_THREADS=1 \
SIM_OUTPUT=outcome/checks/Enterobacteriaceae \
Rscript scripts/simulation/CRC_Entero_Sim.R
```

Use a separate output directory for checks: the workflow updates its setting
and combined output files. To inspect the loaded settings and functions in
RStudio without running the loop, set
options(crc.simulation.source.only = TRUE) before sourcing an entrypoint.

## Outputs

Outputs are created under outcome/simulation/{power,type1}/<family>/:

- *_pvalues.csv: fixed replicate-by-method rows with resampling p-values,
  asymptotic p-values where available, seeds, errors/warnings, sample/cluster
  counts, requested/completed draws and runtime;
- *_summary.csv: rejection rates, 95% Wilson Monte Carlo intervals, failures,
  valid/attempted replicate counts and budget/runtime summaries;
- *_failures.csv: complete method failures and partial study failures for P;
- *_results.rds: results and configuration;
- configuration.rds: settings, fitted GD parameters, input/source checksums,
  package/dependency versions, RNG specification and session information.

Failed p-values remain NA without shifting matrix positions. The primary rate
is the number of valid p-values below 0.05 divided by all planned replicates:
failed replicates are nonrejections. Valid-only rates and failure rates are
also reported. Fisher combination retains the original available-study rule;
the number of valid studies and any partial failure are recorded separately.

n_perm_configured is a requested per-fit budget. n_perm_tests is five for P
and one for CD/RD. Successful DM permutation draws are counted explicitly
(summed over studies for P). The miLineage p-value interfaces do not expose
actual completed draw counts; those entries are NA, not fabricated budgets.
MetaQCAT completed draws are read from its bootstrap arrays.

Elapsed runtimes cannot be summed over methods: generation is shared, FE/RE
share null fits, omnibus summaries reuse component times, and CD tests share
preprocessing. Replicate elapsed times and setting wall times are recorded
separately. Execution checks or finite power values do not establish null
calibration; the null simulation is a separate experiment.

## Power visualization

```sh
Rscript scripts/visualization/simulation/stage_power_summaries.R
Rscript scripts/visualization/simulation/simu_power_bar_plot.R
```

Staging copies the two <family>_power_summary.csv files into res/; overwriting
existing summaries requires --overwrite. The plot reads all 18 methods from
these unified summaries. Its two-family-by-four-setting layout, colors,
type sizes and bottom legend are unchanged. PDF, PNG and plot-data CSV files
are written to plots/. Plotting does not rerun statistical analyses.
