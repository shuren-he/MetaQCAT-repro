# MetaQCAT: paper code and R package

Simulation, CRC/IBD analysis, comparator methods, and figure code for the
MetaQCAT paper. The only project package included is **MetaQCAT 0.2.1**, using
cluster-level wild score bootstrap. Package-bundled example data are retained.
All ten frozen paper inputs are included, including both fitted simulation
GD parameter objects and the CRC/IBD source tables. Generated analysis and
simulation results are not bundled.

## Setup

Use this directory as the working directory in RStudio or the terminal.

```sh
Rscript scripts/check_release.R
Rscript install.R
Rscript scripts/inputs/prepare_inputs.R
```

The preparation command copies and validates the bundled inputs at the
canonical paths without downloading data, fitting models or running analyses.
Public-source reconstruction and provenance are documented in
[docs/INPUTS.md](docs/INPUTS.md).

The installation command installs MetaQCAT, the supplied HMP/miLineage comparator packages, and
analysis/plotting dependencies into `.library/`. R source-package installation
requires a C++ toolchain. Add `--dev` only for package testing/documentation
tools. Installation does not run analyses.

All commands below are run from the repository root. Optional
`rebuild_crc_inputs.R` and `rebuild_ibd_inputs.R` scripts under
`scripts/inputs/` reconstruct source tables from the cited public datasets
into separate audit directories, checking them against the supplied inputs.

## Real-data analyses

CRC and IBD deliberately have **separate entry files**. Each reads its inputs
once and uses one prepared dataset for its MetaQCAT and comparator analyses.
The files are divided into ordinary R sections for interactive execution.

```sh
Rscript scripts/real_data/crc/CRC_all_methods.R
Rscript scripts/real_data/ibd/IBD_all_methods.R
```

| Analysis | Shared preprocessing | Methods |
| --- | --- | --- |
| CRC | Sample depth >=1,000; within-study zero fraction <=0.90; disease-only design | Six MetaQCAT tests; P/CD/RD strategies for QCAT1, QCAT2, DM and ZIGDM |
| IBD | Sample depth >3,000; ten nonoverlapping subject strata; selected-covariate complete cases; zero fraction <=0.90 | Six MetaQCAT tests; P/CD/RD strategies for QCAT.Cluster and QCAT_GEE.Cluster |

Both MetaQCAT analyses use `min.prev = 0.1`, 1,000 bootstrap draws and seed
1234. CRC comparators use a 1,000-permutation budget. IBD comparators use
native between-cluster permutation (`perm.type = "BTW"`) with a requested
cap of 300 permutations per fit. Native adaptive stopping can use fewer draws;
the API does not expose actual completed permutation counts.
CRC is unadjusted; IBD uses the study-specific disease/adjustment design
and subject clusters.
P/CD/RD retain their respective study-level or pooled taxon-screening strategy.
The scripts write fresh p/q-value tables and significant-lineage sets without
requiring previous result-table templates. Running them replaces their outputs.

For IBD, the two comparator tests are `miLineage::QCAT.Cluster` and
`miLineage::QCAT_GEE.Cluster`, labelled QCAT and QCAT2 in the combined table.
P combines study-specific p-values using Fisher's method; CD uses MMUPHin
batch-corrected pooled data and RD uses original pooled data. No DM, ZIGDM
or ConQuR is used in the IBD workflow. Sample selection is shared; P uses
the within-stratum genus filter, whereas CD/RD apply the same zero-fraction
threshold to the pooled original counts. See `docs/INPUTS.md` for details.

By default, the IBD entrypoint refits both MetaQCAT omnibus models (which
provide all six component/omnibus tests) and runs all six comparator analyses.
The complete native computations can take hours. To reuse saved primary
MetaQCAT fits or to validate preprocessing without fitting models:

```sh
IBD_RUN_METAQCAT=false Rscript scripts/real_data/ibd/IBD_all_methods.R
IBD_PREPARE_ONLY=true Rscript scripts/real_data/ibd/IBD_all_methods.R
```

Fit reuse requires matching prepared data, bootstrap budget, seed, package
version and cluster information. Preparation-only mode saves the prepared
inputs without fitting models.

Sensitivity analyses use the prepared primary-analysis objects:

```sh
Rscript scripts/sensitivity/CRC_leave_one_study_out.R
Rscript scripts/sensitivity/IBD_leave_one_cohort_out.R
```

CRC uses 1,000 bootstrap draws. IBD uses saved study-level scores and a matched
300-draw baseline/omission comparison; its primary 1,000-draw results remain
separate.

## Simulations

Simulation inputs are read from `data/simu/`, or automatically from an existing
`simulation/data/` directory when `data/simu/` is absent. No inputs are copied
or deleted automatically.

```sh
Rscript scripts/simulation/CRC_Entero_Sim.R
Rscript scripts/simulation/CRC_Pepto_Sim.R
Rscript scripts/simulation/CRC_Entero_Sim_cluster.R
Rscript scripts/simulation/CRC_Pepto_Sim_cluster.R
Rscript scripts/simulation/CRC_Entero_Type1.R
Rscript scripts/simulation/CRC_Pepto_Type1.R
```

Each replicate generates data once. The original `CRC_*_Sim.R` entrypoints
and independent null retain all 18 methods: six MetaQCAT tests plus P/CD/RD
with QCAT1, QCAT2, DM and ZIGDM. Their independent-sample settings are unchanged.
The new `CRC_*_Sim_cluster.R` entrypoints run correlated power.
Both correlated power and correlated null use 12 methods: the same six
MetaQCAT tests plus P/CD/RD with native `QCAT.Cluster` and
`QCAT_GEE.Cluster`, passing subject IDs explicitly, without DM or ZIGDM.
CD uses MMUPHin only. Correlated power uses the same GD-pair generator and
the same signal mechanisms and strengths as independent power.

Defaults are 1,000 power replicates, 5,000 type-I-error replicates per design,
300 fixed MetaQCAT bootstrap draws and 300 requested permutations per
comparator fit, including RD-DM. Actual adaptive native draw counts are
unavailable and recorded as NA.
Both power designs have four settings per family. Independent-power outputs
remain under `outcome/simulation/power/<family>/`; correlated-power outputs
use the separate `outcome/simulation/power_cluster/<family>/` directory.
Type I error runs both independent and correlated designs by default.
To run only the correlated null design:

```sh
SIM_TYPE1_DESIGNS=correlated Rscript scripts/simulation/CRC_Entero_Type1.R
SIM_TYPE1_DESIGNS=correlated Rscript scripts/simulation/CRC_Pepto_Type1.R
```

`SIM_TYPE1_DESIGNS` accepts `independent`, `correlated`, or a comma-separated
list. Existing `SIM_REPS`, `SIM_BOOTS`, `SIM_PERMS`, `SIM_THREADS` and
`SIM_OUTPUT` controls remain available. Native BTW calls use sorted cluster
rows and study-prefixed pooled IDs: P permutes within study; CD/RD exchange
labels across all pooled clusters without extra study blocking. Per-setting
outputs contain 18 or 12 methods as appropriate, and configuration records
`methods.by.setting`. See
[docs/SIMULATION_WORKFLOWS.md](docs/SIMULATION_WORKFLOWS.md) for budgets,
controls and p-value, Monte Carlo interval, failure and runtime outputs.

## Figures

After the corresponding analyses finish:

```sh
Rscript scripts/visualization/simulation/stage_power_summaries.R
Rscript scripts/visualization/simulation/simu_power_bar_plot.R
Rscript scripts/visualization/simulation/simu_cluster_power_bar_plot.R
Rscript scripts/visualization/crc/CRC_upset_plot.R
Rscript scripts/visualization/crc/CRC_loso_stability_plot_bootstrap.R
Rscript scripts/visualization/ibd/IBD_upset_plot_bootstrap.R
Rscript scripts/visualization/ibd/IBD_loco_stability_plot_bootstrap.R
```

The cluster-power plot reads the two correlated-power summaries directly from
`outcome/simulation/power_cluster/<family>/` and saves separate
`plots/simu_res_cluster_power.*` outputs. It includes all 12 methods;
QCAT1 and QCAT2 denote `QCAT.Cluster` and `QCAT_GEE.Cluster`, respectively.
The independent-power staging and plotting commands above are unchanged.

The IBD UpSet groups are RD, P, CD, FE and RE, as in CRC. Each of RD/P/CD
is the union of discoveries from its two cluster-QCAT tests at permutation
q < 0.05; FE and RE use only their bootstrap omnibus discoveries.

The summary plots retain the paper's layout, including forest panels:

```sh
CRC_PLOT_STYLE=logratio Rscript scripts/visualization/crc/CRC_plot_selected_O_lineages.R
Rscript scripts/visualization/ibd/IBD_extract_selected_O_lineages_latest.R
IBD_PLOT_STYLE=logratio IBD_PLOT_LAYOUT=full Rscript scripts/visualization/ibd/IBD_plot_selected_O_lineages.R f__Peptostreptococcaceae
IBD_PLOT_STYLE=logratio IBD_PLOT_LAYOUT=four_row Rscript scripts/visualization/ibd/IBD_plot_selected_O_lineages.R f__Ruminococcaceae
```

Use `CRC_PLOT_STYLE=violin` or `IBD_PLOT_STYLE=violin` for the alternative
proportion displays. CRC defaults to Enterobacteriaceae, Lachnospiraceae,
Peptostreptococcaceae and Clostridiaceae; family names can also be supplied as
command-line arguments.

## Source checks

Run this on a clean upload tree, before input preparation or analyses:

```sh
Rscript scripts/check_release.R
```

It parses R files, verifies the ten bundled inputs and excludes runtime copies,
results and build artifacts from the release. After preparing inputs, validate
the working analysis inputs with
`Rscript scripts/inputs/prepare_inputs.R --check-only` instead.
These checks do not fit models.
Source checks and preprocessing checks do not replace full
analysis reruns. This organization step does not launch those reruns or upload
anything to GitHub.

See [CODE_INDEX.md](CODE_INDEX.md) for the code map and
[docs/THIRD_PARTY.md](docs/THIRD_PARTY.md) for licenses and attribution.
