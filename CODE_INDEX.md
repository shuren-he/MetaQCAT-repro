# Paper code map

| Location | Purpose |
| --- | --- |
| `packages/MetaQCAT/` | Current MetaQCAT 0.2.1 source, documentation, tests and example data |
| `vendor/HMP/`, `vendor/miLineage/` | Supplied comparator dependency sources |
| `scripts/real_data/crc/CRC_all_methods.R` | Shared CRC preprocessing, all paper methods and comparison tables |
| `scripts/real_data/ibd/IBD_all_methods.R` | Shared IBD preprocessing, all paper methods and comparison tables |
| `scripts/sensitivity/CRC_leave_one_study_out.R` | CRC full-data and leave-one-study-out fits |
| `scripts/sensitivity/IBD_leave_one_cohort_out.R` | IBD cached-score, matched-budget leave-one-cohort-out analysis |
| `scripts/simulation/CRC_Entero_Sim.R`, `CRC_Pepto_Sim.R` | Power entrypoints, all methods on shared generated data |
| `scripts/simulation/CRC_Entero_Type1.R`, `CRC_Pepto_Type1.R` | Independent/correlated null entrypoints, all methods |
| `scripts/simulation/simulation_workflow.R` | Common simulation settings, loops and outputs |
| `scripts/simulation/simulation_helpers.R` | GD generation, signals and method calculations |
| `scripts/visualization/simulation/stage_power_summaries.R` | Collect the two unified family power summaries |
| `scripts/visualization/simulation/simu_power_bar_plot.R` | Paper's two-family-row, four-setting-column power bars |
| `scripts/visualization/crc/` | CRC UpSet, stability and lineage summary figures |
| `scripts/visualization/ibd/` | IBD UpSet, stability, score extraction and lineage figures |
| `code/MetaSub.R` | Current one-step score-summary extraction for descriptive figures |
| `code/IBD/QCAT_C.R` | QCAT-C comparator used by the IBD entrypoint |
| `R/project_paths.R` | Repository paths and local dependency library |
| `install.R` | Dependency and source-package installation; no analyses |
| `scripts/check_release.R` | Static syntax, exclusion and checksum checks |

Only paper workflows and the latest project package are included. External
inputs and generated results are described in `docs/INPUTS.md` and are not
bundled.
