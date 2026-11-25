# Paper code map

| Location | Purpose |
| --- | --- |
| `scripts/inputs/prepare_inputs.R` | Prepare and validate all ten bundled frozen paper inputs; no downloads or fitting |
| `scripts/inputs/rebuild_crc_inputs.R` | Public CRC profiles + supplied frozen taxonomy mapping to verified count/metadata inputs |
| `scripts/inputs/rebuild_ibd_inputs.R` | Frozen public IBD genus-level object to three exactly matching input objects |
| `real/CRC/data/`, `real/IBD/data/`, `simulation/data/` | Ten frozen source input files, including both GD model objects |
| `packages/MetaQCAT/` | Current MetaQCAT 0.2.1 source, documentation, tests and example data |
| `vendor/HMP/`, `vendor/miLineage/` | Supplied comparator dependency sources |
| `scripts/real_data/crc/CRC_all_methods.R` | Shared CRC preprocessing, all paper methods and comparison tables |
| `scripts/real_data/ibd/IBD_all_methods.R` | Shared IBD preprocessing, six MetaQCAT tests, P/CD/RD cluster-QCAT comparisons and tables |
| `scripts/sensitivity/CRC_leave_one_study_out.R` | CRC full-data and leave-one-study-out fits |
| `scripts/sensitivity/IBD_leave_one_cohort_out.R` | IBD cached-score, matched-budget leave-one-cohort-out analysis |
| `scripts/simulation/CRC_Entero_Sim.R`, `CRC_Pepto_Sim.R` | Independent-power entrypoints: original 18 methods on shared generated data |
| `scripts/simulation/CRC_Entero_Sim_cluster.R`, `CRC_Pepto_Sim_cluster.R` | Correlated-power entrypoints: six MetaQCAT and six P/CD/RD cluster-QCAT methods on shared data |
| `scripts/simulation/CRC_Entero_Type1.R`, `CRC_Pepto_Type1.R` | Null entrypoints: independent 18 methods; correlated 12 methods; both designs by default |
| `scripts/simulation/simulation_workflow.R` | Common settings, design selection (`SIM_TYPE1_DESIGNS`), loops and design-specific outputs/configuration |
| `scripts/simulation/simulation_helpers.R` | Unchanged GD/Poisson pair generation, signals, MetaQCAT and design-specific comparator calculations |
| `scripts/visualization/simulation/stage_power_summaries.R` | Collect the two unified family power summaries |
| `scripts/visualization/simulation/simu_power_bar_plot.R` | Paper's two-family-row, four-setting-column power bars |
| `scripts/visualization/simulation/simu_cluster_power_bar_plot.R` | Correlated-power bars read directly from the separate `power_cluster` summaries |
| `scripts/visualization/crc/` | CRC UpSet, stability and lineage summary figures |
| `scripts/visualization/ibd/` | IBD five-group RD/P/CD/FE/RE UpSet, stability, score extraction and lineage figures |
| `code/MetaSub.R` | Current one-step score-summary extraction for descriptive figures |
| `R/project_paths.R` | Repository paths and local dependency library |
| `install.R` | Dependency and source-package installation; no analyses |
| `scripts/check_release.R` | Static syntax and file-exclusion checks |

This index lists the current paper workflows and latest project package.
All required frozen inputs, public-source reconstruction commands and
provenance are described in `docs/INPUTS.md`. Generated results are not bundled.

Independent-null/power comparisons retain P/CD/RD with QCAT1, QCAT2, DM and
ZIGDM. Correlated null and power instead use P/CD/RD with native `QCAT.Cluster` and
`QCAT_GEE.Cluster`, alongside the six MetaQCAT tests; no DM or ZIGDM is run
for either correlated design. CD uses MMUPHin only. See `docs/SIMULATION_WORKFLOWS.md` for
native BTW cluster handling, requested versus completed draws, the correlated-
only command and the `methods.by.setting` output contract.
