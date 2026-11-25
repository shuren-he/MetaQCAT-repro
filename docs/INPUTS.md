# Paper inputs, provenance and preparation

The repository includes all ten frozen input files required by the current
paper workflows: six CRC/IBD source tables and four simulation inputs,
including both fitted generalized-Dirichlet (GD) parameter objects. They are
inputs, not fitted association results or simulation outputs. Package examples
are included separately and are not substitutes for these full inputs.

Run this once from the repository root, before starting analyses:

```sh
Rscript scripts/inputs/prepare_inputs.R
Rscript scripts/inputs/prepare_inputs.R --check-only
```

The first command copies the supplied files from `real/CRC/data/`,
`real/IBD/data/` and `simulation/data/` into `data/CRC/`, `data/IBD/` and
`data/simu/`. It validates object dimensions, sample/taxon alignment, count
values and GD parameter structure. It uses base R, makes no network requests,
and does not fit models. Repeating it is safe; different existing files are
not overwritten unless `--overwrite` is explicitly supplied. Canonical
`data/` files are runtime copies ignored by Git, not missing release assets.

## CRC

These frozen files are supplied under `real/CRC/data/`; the preparation
command copies them to `data/CRC/`:

| File | Required object |
| --- | --- |
| `count.genus.RData` | `count.genus`: sample-by-genus integer count matrix |
| `meta.RData` | `meta`: metadata containing `Sample_ID`, `Country` and `Group` (`CRC`/`CTR`) |
| `tax.RData` | `tax`: taxonomy with kingdom through genus columns |

The count row names must match `make.names(meta$Sample_ID)`. Genus columns
must match the taxonomy's genus identifiers. CRC preprocessing retains 573
samples in five studies for the paper inputs; the depth rule removes one
sample, not an entire study.

The source is Wirbel et al. (2019), fixed Zenodo deposit
[3517209, version 1.0](https://doi.org/10.5281/zenodo.3517209), with
`species_profiles_g2_l75_motus2.0.0.tsv` and `meta_all.tsv`.
The supplied 574-by-133 integer matrix is an exact genus aggregation of
count-scale mOTUs marker profiles, not relative profiles multiplied by an
assumed library size. Its depths are reported marker-profile count depths.

An optional public-source reconstruction is provided:

```sh
Rscript scripts/inputs/rebuild_crc_inputs.R
```

This needs `data.table` (installed by `install.R`). It downloads the two fixed
source files, or reuses a verified cache, and writes only to
`data/rebuilt/CRC/`. The feature screen is applied to the original 575
CRC/control samples from AT-CRC, CN-CRC, FR-CRC, DE-CRC and US-CRC:
a maximum relative abundance of at least 0.001 in at least three studies,
excluding the unassigned feature. The 849 selected mOTUs are summed on their
original count scale using the supplied `tax.RData` genus mapping.

The supplied taxonomy is a frozen project annotation mapping; it is distributed
as an input rather than claimed to be inferred from the count matrix. The
reconstruction checks its mOTU selection against the public source screen.
It also reproduces all 18 metadata columns, including `block`.
Counts are identical to the frozen matrix; all metadata column values and
types match (the rebuilt metadata uses a plain data.frame).

The frozen 574-sample universe excludes `ERR1018294` from the 575 source
samples. The historical reason is not recorded; the reconstruction preserves
this explicit choice and does not mislabel it as a depth exclusion. The later
analysis depth rule removes `CCIS12370844ST-4-0` (depth 611), leaving 573.
AUS is the existing code for the Austrian AT-CRC study, not Australia.
These preparation scripts do not change either sample-selection rule.

## IBD

These frozen files are supplied under `real/IBD/data/`; the preparation
command copies them to `data/IBD/`. The IBD entrypoint also accepts the
supplied location directly if the canonical directory is absent:

| File | Required object |
| --- | --- |
| `count.genus.RData` | `otu`: genus-by-sample count matrix |
| `meta.RData` | `meta.data`: metadata aligned to `colnames(otu)` by row name |
| `tax.RData` | `Tax`: aligned genus taxonomy, first six ranks |

The source is Ma et al. (2022), [Population structure discovery in
meta-analyzed microbial communities and inflammatory bowel disease](https://doi.org/10.1186/s13059-022-02753-4).
The authors' [IBD analysis repository](https://github.com/biobakery/ibd_paper)
provides `data/physeq/genera.RData`. We fix the source to commit
`ae9d98ad170546f02af3c67797cdba971bc95e7d`; the exact downloadable URL
is in the reconstruction script:

```sh
Rscript scripts/inputs/rebuild_ibd_inputs.R
```

It downloads the publication's already-cleaned, pre-filtered genus-level
phyloseq object and extracts its count, sample-data and taxonomy slots.
No extra filtering, normalization, imputation or reordering is applied.
The resulting `otu` (249-by-4,789), `meta.data` (4,789-by-57) and
`Tax` (249-by-8) are each identical to the corresponding frozen object,
including types and dimnames. No phyloseq installation is needed for this
fixed serialized source. Earlier sequencing/OTU processing is upstream work,
documented in that repository's `rmds/` scripts; this route starts at its
published genus-level count object. Rebuilt files go to `data/rebuilt/IBD/`,
not over the primary analysis inputs.

Metadata fields used are `dataset_name`, `subject_accession`, `sample_type`,
`disease`, `gender`, `race` and `age`. The primary workflow produces ten
nonoverlapping strata from seven source cohorts: 3,141 samples and 1,513
subjects for the paper inputs. Sample exclusions and study summaries are
exported by `IBD_all_methods.R`.
The older `real/IBD/code/IBD_process.R` entry forwards to this same workflow
when run from the repository root.

The primary sample rule is raw depth strictly greater than 3,000.
Complete cases are required for the adjustment covariates selected in each
stratum; participants are retained in only one stratum. MetaQCAT and P use
genera with within-stratum zero fraction <=0.90. CD/RD use the same retained
samples but screen the pooled original count matrix at zero fraction <=0.90,
rather than zero-filling taxa already removed by the study-specific screen.

The IBD comparators are `QCAT.Cluster` and `QCAT_GEE.Cluster` under each of
P, CD and RD (combined-table labels QCAT and QCAT2). Disease is the tested variable;
the selected stratum-specific nuisance covariates are retained, without
adding study-indicator columns. P uses each study's adjusted design.
Pooled cluster IDs are prefixed by study to keep different subjects distinct.
The two-part QCAT2 input is rarefied within each study for P, or on the pooled
matrix for CD/RD. CD performs MMUPHin correction only. Both native tests
use `perm.type = "BTW"` and a default requested cap of 300 permutations per
fit; their actual adaptive permutation counts are not exposed by the API.

The IBD entrypoint saves comparator checkpoints, p/q-values and call diagnostics
under `res/IBD/cluster_qcat/`. The updated
`tables/IBD_family_q_val.csv` and
`tables/IBD_family_q_val_resampling.csv` include the current MetaQCAT and
cluster-QCAT results; the resampling table has six MetaQCAT bootstrap and
six comparator permutation columns. The canonical
`res/IBD/IBD_sig_family.csv` contains RD, P, CD, FE and RE.
RD/P/CD mark the union of their two comparator discoveries at q < 0.05;
FE/RE mark only FE-MetaQCAT-O and RE-MetaQCAT-O discoveries.
The UpSet script writes `outcome/IBD/IBD_overlap_plot.pdf` and its PNG
preview without rerunning statistical tests.

`res/IBD/IBD_lineage_combine.csv` now has 17 named columns: lineage,
the FE/RE omnibus p/q-values, and p/q-values for the six comparator tests.
Summary-plot extraction reads the primary FE/RE RDS objects directly; it
does not require the former 15-column QCAT-C combined table.

By default, `IBD_all_methods.R` runs primary MetaQCAT and comparator fits.
`IBD_RUN_METAQCAT=false` reuses saved primary MetaQCAT fits only after
checking prepared-data equality, bootstrap budget, seed, package version
and cluster information. `IBD_PREPARE_ONLY=true` saves the prepared inputs
and exclusion reports without fitting models.
These controls do not download inputs or regenerate missing fitted results;
use the input-preparation command above for the supplied source tables.

The reference IBD run used R 4.6.0, MetaQCAT 0.2.1, miLineage 4.0, MMUPHin 2.0.0
and digest 0.6.39. The entrypoint fixes the Mersenne-Twister RNG with
Inversion/Rejection sampling and records input/script MD5 values, package
versions and session information alongside the outputs.

## Simulation

All four frozen inputs are supplied under `simulation/data/`, including the
two small fitted GD objects. The preparation command copies them to
`data/simu/`; simulation entrypoints also accept the supplied location
directly if the canonical directory is absent:

- `mod_gdm_lst_Enterobacteriaceae.RData` and
  `mod_gdm_lst_Peptostreptococcaceae.RData`, each containing `mod_gdm_lst`;
- `meta.rda`, containing `meta`;
- `count.genus.rda`, containing sample-by-genus `count.genus`.

The simulation count matrix and metadata are identical to the supplied CRC
source objects, before the analysis-level depth filter. They therefore have
the same public source and preparation route described above.

Each GD file contains `mod_gdm_lst`, a named list ordered AUS, CHI, FRA,
GER, USA. Each study has `est.para` (including `alpha.est`, `phi.est`,
convergence status and iterations) and `col.order`. The Enterobacteriaceae
file is 695 bytes: four genera, with a 2-by-3 alpha matrix per study.
The Peptostreptococcaceae file is 548 bytes: three genera, with a 2-by-2
alpha matrix per study. Taxon ordering and exact genus identifiers are fixed
in the family entrypoints; all saved column-order mappings are identity.

These are the original project parameter objects used by the paper workflow,
not new fits made during release preparation. Historical project fitting code
used study-specific `GDM.reg` models with intercept and disease columns,
fitting beta-binomial sticks through an EM/BFGS procedure. The files do not
store their exact fit calls or software session, so we do not claim that a
fresh optimization would reproduce every saved parameter. The complete frozen
objects are supplied instead; reproducing the advertised simulations does
not require refitting them.

Current baseline generation evaluates the saved model at `X = (1, 0)`,
using the intercept row of `alpha.est` and the saved dispersion parameters.
Power signals are added afterward. The stored disease coefficients do not
introduce a baseline case/control signal. See `SIMULATION_WORKFLOWS.md`
for the independent/correlated generators, settings and output schemas.

## Optional source reconstruction controls

Both public-source scripts support `INPUT_SOURCE_CACHE` for existing source
downloads and `INPUT_REBUILD_DIR` for a separate reconstruction directory.
They verify reconstructed objects against the supplied frozen references.
The CRC taxonomy remains the supplied frozen annotation mapping; a new
annotation database is not substituted. CRC reconstruction requires
`data.table`; IBD reconstruction uses base R.

The public-source reconstructions are optional audits, not prerequisites
for the normal frozen-input route. Their outputs do not replace `data/CRC/`
or `data/IBD/` automatically, and neither script refits GD parameters or
runs association tests.

## Outputs and figure dependencies

Real-data entrypoints create prepared inputs, fits and tables under `res/`
and `tables/`. The sensitivity and visualization scripts consume those files;
they do not require manually edited old p/q-value tables. Generated figures
are written under `plots/` or `outcome/`. Simulation outputs are described in
`SIMULATION_WORKFLOWS.md`.

Full-analysis input files and package examples are bundled. Generated result
objects, public-source download caches and PDFs are not bundled.
The exact source assembly of
the supplementary source-composition figure (`composition_info.pdf`) was
not found, so no substitute
script is presented as an exact reproduction of that figure.
