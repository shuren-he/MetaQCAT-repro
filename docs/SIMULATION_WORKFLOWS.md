# Paper simulation workflows

Run from the repository root after installing the bundled MetaQCAT 0.2.1 and
comparator dependencies (HMP, miLineage and MMUPHin). There is one generation
and preprocessing workflow, not separate generators for our tests and baselines.

## Bundled frozen inputs

All four required files are included under `simulation/data/`.
Prepare their canonical paths once from the repository root:

```sh
Rscript scripts/inputs/prepare_inputs.R
```

The simulation scripts read the resulting `data/simu/` files; if that
directory is absent, they also accept the bundled `simulation/data/` files:

- mod_gdm_lst_Enterobacteriaceae.RData
- mod_gdm_lst_Peptostreptococcaceae.RData
- meta.rda, containing meta
- count.genus.rda, containing the aligned count.genus table

The model files contain the five study-specific fitted generalized Dirichlet
parameter sets and taxon-order mappings used by the frozen experiment.
They are supplied as inputs, so no fitting script or separate data request
is needed to run the simulations. The simulation count/metadata objects
are identical to the frozen CRC source tables. Public source links,
source-to-input reconstruction scripts and model provenance are documented
in [INPUTS.md](INPUTS.md). The release preserves the saved GD parameters;
it does not silently substitute a new fit.

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
Rscript scripts/simulation/CRC_Entero_Sim_cluster.R
Rscript scripts/simulation/CRC_Pepto_Sim_cluster.R
Rscript scripts/simulation/CRC_Entero_Type1.R
Rscript scripts/simulation/CRC_Pepto_Type1.R
```

Each replicate generates one dataset for every method in its setting.
Independent-sample power and independent type I error retain the original
18 methods and their existing settings:

- FE/RE-MetaQCAT, FE/RE-MetaQCAT* and FE/RE-MetaQCAT-O;
- P, CD and RD, each with QCAT1, QCAT2, DM and ZIGDM.

Correlated-sample power and correlated type I error instead evaluate 12 methods:

- the same six MetaQCAT component/omnibus tests;
- P, CD and RD, each with `miLineage::QCAT.Cluster` and
  `miLineage::QCAT_GEE.Cluster`.

DM and ZIGDM are not run in either correlated design. The two new
`CRC_*_Sim_cluster.R` entrypoints run correlated power only; the existing
`CRC_*_Sim.R` entrypoints still run independent power. Their inputs, signals,
budgets, method set and output locations are unchanged.

Every method receives the same generated counts, labels, retained sample rows
and taxon identities. P combines single-study p-values using Fisher's method.
CD pools data after MMUPHin batch adjustment only, preserving the simulated
case/control variable; no ConQuR correction is used. RD pools the generated
data directly. Method-specific rarefaction and internal eligibility rules
remain part of each test: QCAT2 and `QCAT_GEE.Cluster` use once-rarefied data
(separately per study for P, or per pooled CD/RD table). The other tests retain
their corresponding strategy inputs, rather than all receiving rarefied data.

The six entrypoints set family constants and target, then source
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
| Independent-null/power comparator budget per fit | 300 for all tests, including RD-DM |
| Correlated-sample native permutation cap per fit | 300; native adaptive stopping |
| Requested workers | 7; capped by available cores, 1 on Windows |
| Significance level | 0.05 |
| MetaQCAT minimum prevalence | 0.1 |

Both independent and correlated power have two mechanisms (proportional and
count-addition) crossed with two patterns (homogeneous and one-study-reversed),
with 1,000 replicates per setting and family. One taxon is randomly chosen.
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
signal size zero. `SIM_TYPE1_DESIGNS` can select either design; it does not
change the generator. Correlated samples are randomly paired within study. Both
members' independent GD vectors are averaged with a third, shared GD baseline
vector. Labels are assigned to clusters. An odd remaining sample stays in a
singleton cluster. MetaQCAT uses these subject IDs for the wild bootstrap.
The study-mean Poisson depths and GD/shared-vector pairing mechanism are
unchanged by the comparator update.

Correlated power uses the same shared-GD pairing and cluster-level case/control
assignment. Signals are added using the independent-power procedure and
strengths above: proportional signals modify the averaged proportions before
multinomial sampling; count signals are added after multinomial sampling and
rescaled to each sample's original simulated depth. There is no second
data-generation call for a different method.

Native QCAT calls in both correlated designs use between-cluster permutation
(`perm.type = "BTW"`). Input rows are sorted into contiguous subject clusters,
and pooled subject IDs are prefixed by study. P fits each study separately, so
its cluster-label permutations stay within study. CD/RD exchange labels across
all pooled clusters, with no additional study blocking. Unequal cluster sizes
are retained; neither individual rows nor cluster sizes are shuffled.
Both `QCAT.Cluster` and `QCAT_GEE.Cluster` receive the matched subject vector
explicitly as `ID = id`, with `perm.type = "BTW"`. The calls to the older
`QCAT` and `QCAT_GEE` interfaces in `simulation_helpers.R` are retained only
for the independent-sample comparator branch, not for correlated simulations.

For MetaQCAT, each count/zero null model is fitted once per study and replicate.
FE/RE reuse the scores and the count/zero components share keyed cluster
weights. Marginal p-values are (1 + exceedances)/(n.boot + 1); -O uses the
current package's Cauchy combination of these marginal p-values. RE count and
zero components both use the positive part of the summed variance-component
score. DM, used only for independent-null/power comparisons, preserves its
exceedances/successful-permutations rule without a plus-one correction.
miLineage retains its own permutation backend, including adaptive stopping
for the native cluster tests. MetaQCAT uses all 300 requested bootstrap draws
by default, without adaptive stopping. The configured budgets are therefore
not a single common resampling scheme.

## Small checks and controls

Environment controls are SIM_REPS, SIM_BOOTS, SIM_PERMS, SIM_RD_DM_PERMS,
SIM_THREADS and SIM_OUTPUT. They apply to power and type I error entrypoints.
SIM_RD_DM_PERMS affects only the independent-null/power RD-DM comparison
and defaults to 300, matching the other comparator budgets.
For type I error, SIM_TYPE1_DESIGNS accepts `independent`, `correlated`, or the
comma-separated list `independent,correlated` (the default).

```sh
SIM_REPS=2 SIM_BOOTS=19 SIM_PERMS=19 SIM_RD_DM_PERMS=19 SIM_THREADS=1 \
SIM_OUTPUT=outcome/checks/Enterobacteriaceae \
Rscript scripts/simulation/CRC_Entero_Sim.R
```

Run only the correlated null design for either family with:

```sh
SIM_TYPE1_DESIGNS=correlated Rscript scripts/simulation/CRC_Entero_Type1.R
SIM_TYPE1_DESIGNS=correlated Rscript scripts/simulation/CRC_Pepto_Type1.R
```

Run a small correlated-power check without running independent power:

```sh
SIM_REPS=2 SIM_BOOTS=19 SIM_PERMS=19 SIM_THREADS=1 \
SIM_OUTPUT=outcome/checks/Enterobacteriaceae_correlated_power \
Rscript scripts/simulation/CRC_Entero_Sim_cluster.R
```

`CRC_Pepto_Sim_cluster.R` uses the same controls for the other family.
`SIM_TYPE1_DESIGNS` affects type I error only; it is not needed for the
correlated-power entrypoints.

Use a separate output directory for checks: the workflow updates its setting
and combined output files. To inspect the loaded settings and functions in
RStudio without running the loop, set
options(crc.simulation.source.only = TRUE) before sourcing an entrypoint.

## Outputs

Outputs are created under `outcome/simulation/power/<family>/` for independent
power, `outcome/simulation/power_cluster/<family>/` for correlated power, and
`outcome/simulation/type1/<family>/` for type I error. Correlated-power runs
therefore do not overwrite independent-power results.

The correlated-power combined files are
`<family>_correlated_power_resampling_pvalues.csv`,
`<family>_correlated_power_resampling_summary.csv`, and
`<family>_correlated_power_results.rds`. Their rows identify
`design = "correlated"`; summary `power` is the same all-attempted
rejection rate, and `n_perm_requested` records the native permutation cap.
These aliases do not change the independent output schema.

- *_pvalues.csv: fixed replicate-by-method rows with resampling p-values,
  asymptotic p-values where available, seeds, errors/warnings, sample/cluster
  counts, requested/completed draws and runtime;
- *_summary.csv: rejection rates, 95% Wilson Monte Carlo intervals, failures,
  valid/attempted replicate counts and budget/runtime summaries;
- *_failures.csv: complete method failures and partial study failures for P;
- *_results.rds: results and configuration;
- configuration.rds: settings, fitted GD parameters, input/source checksums,
  package/dependency versions, RNG specification and session information;
  `methods.by.setting` records the method names used for each setting.

Per-setting tables contain 18 methods for independent power/null and 12 for
correlated power/null. Combined type I error outputs retain these design-specific
method sets; an absent method is not reported as an executed failed test.

Failed p-values remain NA without shifting matrix positions. The primary rate
is the number of valid p-values below 0.05 divided by all planned replicates:
failed replicates are nonrejections. Valid-only rates and failure rates are
also reported. Fisher combination retains the original available-study rule;
the number of valid studies and any partial failure are recorded separately.

n_perm_configured is a requested per-fit budget. n_perm_tests is five for P
and one for CD/RD. Successful DM permutation draws are counted explicitly
(summed over studies for P), only where DM is part of the method set. The
miLineage p-value interfaces do not expose actual completed draw counts;
those entries are NA, not fabricated budgets. In particular, 300 is the
requested native cluster-permutation cap, not an asserted completed count.
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
Rscript scripts/visualization/simulation/simu_cluster_power_bar_plot.R
```

Staging copies the two <family>_power_summary.csv files into res/; overwriting
existing summaries requires --overwrite. The plot reads all 18 methods from
these unified summaries. Its two-family-by-four-setting layout and bottom
legend are retained. Italic family labels and two-level Low/High heterogeneity
and Small/Large ZP change headers are included in the figure. Strategy labels
and y-axis ticks are 35% larger; MetaQCAT-O uses deep brown (`#5C3317`).
PDF, PNG and plot-data CSV files
are written to plots/. Plotting does not rerun statistical analyses.

The correlated-power plot reads the two combined summary files directly from
`outcome/simulation/power_cluster/<family>/`; it does not use or replace the
independent staging files. It includes all 12 correlated-sample methods, using
MetaQCAT bootstrap power and cluster-QCAT permutation power. In its legend,
QCAT1 denotes `QCAT.Cluster` and QCAT2 denotes `QCAT_GEE.Cluster`.
Its matching two-by-four layout is saved as `plots/simu_res_cluster_power.pdf`,
`plots/simu_res_cluster_power.png`, and
`plots/simu_res_cluster_power_plot_data.csv`.
