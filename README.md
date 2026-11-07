# MetaQCAT: paper code and R package

Simulation, CRC/IBD analysis, comparator methods, and figure code for the
MetaQCAT paper. The only project package included is **MetaQCAT 0.2.1**, using
cluster-level wild score bootstrap. Package-bundled example data are retained;
external analysis data and generated results are not included.

## Setup

Use this directory as the working directory in RStudio or the terminal.

```sh
Rscript install.R
```

This installs MetaQCAT, the supplied HMP/miLineage comparator packages, and
analysis/plotting dependencies into `.library/`. R source-package installation
requires a C++ toolchain. Add `--dev` only for package testing/documentation
tools. Installation does not run analyses.

Restore the inputs listed in [docs/INPUTS.md](docs/INPUTS.md) before running
analyses. All commands below are run from the repository root.

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
| IBD | Sample depth >3,000; nonoverlapping subject strata; selected-covariate complete cases; zero fraction <=0.90 | Six MetaQCAT tests and QCAT-C |

Both MetaQCAT analyses use `min.prev = 0.1`, 1,000 bootstrap draws and seed
1234. The real-data comparator budget is 1,000 permutations. CRC is unadjusted;
IBD uses the study-specific disease/adjustment design and subject clusters.
P/CD/RD retain their respective study-level or pooled taxon-screening strategy.
The scripts write fresh p/q-value tables and significant-lineage sets without
requiring previous result-table templates. Running them replaces their outputs.

Sensitivity analyses use the prepared primary-analysis objects:

```sh
Rscript scripts/sensitivity/CRC_leave_one_study_out.R
Rscript scripts/sensitivity/IBD_leave_one_cohort_out.R
```

CRC uses 1,000 bootstrap draws. IBD uses saved study-level scores and a matched
300-draw baseline/omission comparison; its primary 1,000-draw results remain
separate.

## Simulations

```sh
Rscript scripts/simulation/CRC_Entero_Sim.R
Rscript scripts/simulation/CRC_Pepto_Sim.R
Rscript scripts/simulation/CRC_Entero_Type1.R
Rscript scripts/simulation/CRC_Pepto_Type1.R
```

Each replicate generates data once for all six MetaQCAT and twelve P/CD/RD
tests. Defaults are 1,000 power replicates, 5,000 type-I-error replicates and
300 MetaQCAT bootstrap draws. Power has four settings; null runs include
independent and correlated samples. Comparator permutation budgets and the
generated p-value, Monte Carlo interval, failure and runtime summaries are
documented in [docs/SIMULATION_WORKFLOWS.md](docs/SIMULATION_WORKFLOWS.md).

## Figures

After the corresponding analyses finish:

```sh
Rscript scripts/visualization/simulation/stage_power_summaries.R
Rscript scripts/visualization/simulation/simu_power_bar_plot.R
Rscript scripts/visualization/crc/CRC_upset_plot.R
Rscript scripts/visualization/crc/CRC_loso_stability_plot_bootstrap.R
Rscript scripts/visualization/ibd/IBD_upset_plot_bootstrap.R
Rscript scripts/visualization/ibd/IBD_loco_stability_plot_bootstrap.R
```

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

```sh
Rscript scripts/check_release.R
```

This parses R files, checks file exclusions and verifies `SHA256SUMS` without
fitting models. Source checks and preprocessing checks do not replace full
analysis reruns. This organization step does not launch those reruns or upload
anything to GitHub.

See [CODE_INDEX.md](CODE_INDEX.md) for the code map and
[docs/THIRD_PARTY.md](docs/THIRD_PARTY.md) for licenses and attribution.
