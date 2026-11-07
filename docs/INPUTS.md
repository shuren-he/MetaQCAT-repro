# External inputs and generated files

Package examples under `packages/MetaQCAT/data/` are included. The following
full-analysis inputs are separate and are not replaced by those examples.

## CRC

Place these files under `data/CRC/`:

| File | Required object |
| --- | --- |
| `count.genus.RData` | `count.genus`: sample-by-genus integer count matrix |
| `meta.RData` | `meta`: metadata containing `Sample_ID`, `Country` and `Group` (`CRC`/`CTR`) |
| `tax.RData` | `tax`: taxonomy with kingdom through genus columns |

The count row names must match `make.names(meta$Sample_ID)`. Genus columns
must match the taxonomy's genus identifiers. CRC preprocessing retains 573
samples in five studies for the paper inputs; the depth rule removes one
sample, not an entire study.

## IBD

Place these files under `data/IBD/`:

| File | Required object |
| --- | --- |
| `count.genus.RData` | `otu`: genus-by-sample count matrix |
| `meta.RData` | `meta.data`: metadata aligned to `colnames(otu)` by row name |
| `tax.RData` | `Tax`: aligned genus taxonomy, first six ranks |

Metadata fields used are `dataset_name`, `subject_accession`, `sample_type`,
`disease`, `gender`, `race` and `age`. The primary workflow produces ten
nonoverlapping strata from seven source cohorts: 3,141 samples and 1,513
subjects for the paper inputs. Sample exclusions and study summaries are
exported by `IBD_all_methods.R`.

## Simulation

Place the frozen GD models and source compositions under `data/simu/`:

- `mod_gdm_lst_Enterobacteriaceae.RData` and
  `mod_gdm_lst_Peptostreptococcaceae.RData`, each containing `mod_gdm_lst`;
- `meta.rda`, containing `meta`;
- `count.genus.rda`, containing sample-by-genus `count.genus`.

See `SIMULATION_WORKFLOWS.md` for generation settings and output schemas.
The fitted GD input objects are external inputs, not bundled results.

## Outputs and figure dependencies

Real-data entrypoints create prepared inputs, fits and tables under `res/`
and `tables/`. The sensitivity and visualization scripts consume those files;
they do not require manually edited old p/q-value tables. Generated figures
are written under `plots/` or `outcome/`. Simulation outputs are described in
`SIMULATION_WORKFLOWS.md`.

No full-analysis data, result objects, caches, or PDFs are distributed in this
code-only directory. The exact source assembly of the supplementary source-
composition figure (`composition_info.pdf`) was not found, so no substitute
script is presented as an exact reproduction of that figure.
