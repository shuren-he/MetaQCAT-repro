# MetaQCAT Reproducibility Code

R code and processed data for the MetaQCAT simulation studies and real-data analyses of colorectal cancer (CRC) and inflammatory bowel disease (IBD).

The repository includes the MetaQCAT source package, fitted simulation models, processed count tables, metadata, and taxonomy tables. Plotting scripts produce **PDF files only**. Generated figures are not bundled with the repository.

## Repository structure

```text
package/
  MetaQCAT_0.0.3.0000.tar.gz
simulation/
  code/
    CRC_Pepto_Sim.R
    CRC_Entero_Sim.R
    bar_plot.R
  data/
    count.genus.rda
    meta.rda
    mod_gdm_lst_Peptostreptococcaceae.RData
    mod_gdm_lst_Enterobacteriaceae.RData
real/
  CRC/
    code/
      CRC_process.R
      composition_plot.R
    data/
      count.genus.RData
      meta.RData
      tax.RData
  IBD/
    code/
      IBD_process.R
    data/
      count.genus.RData
      meta.RData
      tax.RData
  code/
    upset_plot.R
```

The simulation input tables are derived from the CRC data and are used with the fitted models to generate simulated counts. The composition plot uses the observed data in `real/CRC/data/`.

All examples below assume that the working directory is the repository root. Scripts also locate the repository from their own file paths when run with `Rscript` or `source()`.

### Generated output directories

`plot/` is created locally when a plotting script saves a figure; it does not need to exist after cloning the repository. All three plotting scripts create missing output directories recursively before saving their PDFs:

| Plotting script | Generated output location |
| --- | --- |
| `simulation/code/bar_plot.R` | `plot/simulation/*.pdf` |
| `real/CRC/code/composition_plot.R` | `plot/real/CRC_composition.pdf` |
| `real/code/upset_plot.R` | `plot/real/*_overlap.pdf`, or the PDF path supplied to the function |

No manual directory creation is required. Simulation bar plots require simulation result CSVs, and UpSet plots require discovery results. If those inputs are absent, the scripts do not generate the corresponding figures. Generated files under `plot/real/` are excluded from Git by `.gitignore`.

## Installation

Use R with the dependencies below. A C++ compiler is required to install the MetaQCAT source package. Linux or macOS is suitable for parallel simulations; on Windows, set `nthreads <- 1` because the simulation scripts use `parallel::mclapply()`.

### CRAN dependencies

Run in R from the repository root:

```r
install.packages(c(
  "brglm2", "CompQuadForm", "data.table", "dplyr", "MASS",
  "parallelly", "psych", "Rcpp", "RcppArmadillo",
  "doParallel", "foreach", "dirmult", "stringr", "zoo",
  "vegan", "gplots", "rpart", "rpart.plot",
  "ggplot2", "ggh4x", "patchwork", "UpSetR"
), repos = "https://cloud.r-project.org")
```

HMP is loaded by the real-data analysis scripts; it is not needed for simulation. Install HMP 2.0 from the [CRAN archive](https://cran.r-project.org/src/contrib/Archive/HMP/):

```r
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/HMP/HMP_2.0.tar.gz",
  repos = NULL, type = "source"
)
```

### Bioconductor dependencies

Install qvalue for the real-data analysis scripts through [BiocManager](https://www.bioconductor.org/install/).

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install("qvalue")
```

### MetaQCAT source package

Install the bundled MetaQCAT package:

```r
install.packages("package/MetaQCAT_0.0.3.0000.tar.gz",
                 repos = NULL, type = "source")
```

The simulation scripts require only MetaQCAT and its package dependencies. They do not run comparison methods or require miLineage, MMUPHin, or HMP.

For composition plotting alone, only `ggplot2`, `ggh4x`, and `patchwork` are needed. UpSet plotting requires `UpSetR`, and simulation bar plotting requires `ggplot2`.

## Simulation studies

The two scripts generate simulations for Peptostreptococcaceae and Enterobacteriaceae across five studies: AUS, CHI, FRA, GER, and USA. The simulations evaluate only FE-MetaQCAT, FE-MetaQCAT*, FE-MetaQCAT-O, RE-MetaQCAT, RE-MetaQCAT*, and RE-MetaQCAT-O. The data-generation settings and simulation scenarios are retained.

Before running, adjust the settings near the top of each simulation script to suit the available hardware:

| Setting | Default | Meaning |
| --- | --- | --- |
| `seed` | `2333` | Simulation seed |
| `nthreads` | `125` | Number of parallel replicate workers; reduce for a local machine |
| `n.perm` | `300` | Permutation count used by MetaQCAT |
| `simulation.target` | `"power"` | Use `"power"` or `"type1"` |
| `sim.time` | `1000` / `5000` | Replicates for power / Type I error |

Power simulations include homogeneous, two-null-study, and one-flipped-study scenarios, with positive-part (`proportion`) and zero-part (`count`) signals. Type I error simulations use a zero signal in the homogeneous scenario.

Run from a terminal:

```bash
Rscript simulation/code/CRC_Pepto_Sim.R
Rscript simulation/code/CRC_Entero_Sim.R
```

Replicate-level p-value CSVs are written to `simulation/outcome/pos_part_res/`. Despite its name, this directory contains both positive-part and zero-part simulation results. Each asymptotic result CSV contains the three FE-MetaQCAT columns; each permutation result CSV contains the three FE and three RE columns. Separate failure logs are written when a MetaQCAT replicate fails. Running the same settings again overwrites the matching result files.

### Simulation bar plots

After the simulation CSVs have been generated:

```bash
Rscript simulation/code/bar_plot.R
```

The script reads `simulation/outcome/pos_part_res/` and saves PDFs under `plot/simulation/`, for example:

```text
power_proportion_perm_bar_plot.pdf
power_count_perm_bar_plot.pdf
type1_proportion_perm_bar_plot.pdf
```

Only figures with corresponding input results are generated. The default is permutation results (`inference <- "perm"`) with rejection cutoff `alpha <- 0.05`; change `inference` to `"asy"` for asymptotic results. The bars show the three MetaQCAT variants grouped under FE and RE (FE only for asymptotic results). Older CSVs containing comparison methods can also be read: only MetaQCAT columns are selected. Rates are calculated over valid replicates, excluding missing p-values. The plotting script writes no PNGs or summary CSVs.

## Real-data analyses

Run the CRC and IBD analyses in separate R sessions. These scripts begin by clearing the current workspace and retain analysis results in memory rather than automatically exporting result files.

```r
# CRC analysis
source("real/CRC/code/CRC_process.R")
```

```r
# IBD analysis, in a separate session
source("real/IBD/code/IBD_process.R")
```

Both analyses run FE-MetaQCAT-O and RE-MetaQCAT-O with 1,000 permutations. Set `n.threads` in their MetaQCAT calls to suit your machine. CRC results are stored in `CRC_FE` and `CRC_RE`; IBD results are stored in `IBD_FE` and `IBD_RE`. Their corresponding `*_pval` vectors contain lineage p-values. The `*_qval` vectors are calculated using `qvalue::qvalue()` separately for FE and RE and are used to select discoveries for the UpSet plots.

The IBD workflow uses subject IDs to account for repeated measurements, splits HMP2 into pediatric and adult cohorts, and performs between-subject permutations.

### CRC composition plot

This plot can be generated directly from the processed CRC data without running a MetaQCAT analysis:

```bash
Rscript real/CRC/code/composition_plot.R
```

Output: `plot/real/CRC_composition.pdf`.

The single-page PDF contains two panels: Peptostreptococcaceae and Enterobacteriaceae, each faceted by country and CRC/CTR group. Each bar represents one sample, normalized within its family. Samples with zero total family abundance are omitted. The shallow sample excluded in the original workflow is removed by sample ID. Observed genus counts are retained across countries.

### UpSet discovery plots

The paper's UpSet plots compare discovered lineages with **q-value < 0.05**, using q-values calculated by the R `qvalue` package. After running the CRC analysis in R:

```r
source("real/code/upset_plot.R")
res <- list(
  FE = names(CRC_FE_qval)[!is.na(CRC_FE_qval) & CRC_FE_qval < 0.05],
  RE = names(CRC_RE_qval)[!is.na(CRC_RE_qval) & CRC_RE_qval < 0.05]
)
plot_discovery_upset(res, "plot/real/CRC_overlap.pdf")
```

For IBD, run the IBD analysis in a separate session and replace the q-value vectors with `IBD_FE_qval` and `IBD_RE_qval`, saving to `plot/real/IBD_overlap.pdf`.

Each q-value vector must be named by lineage ID. The example selects discoveries directly from the calculated q-values; no additional BH adjustment is applied. You may also supply a named list of discovered lineage IDs directly to `plot_discovery_upset()` to compare additional methods. At least two methods must have discoveries to draw the plot.

An existing binary discovery table can be read instead: the first column must contain lineage IDs, and the remaining named method columns must contain 0/1 indicators. For a table saved at `real/CRC/outcome/sig_family.csv`:

```bash
Rscript real/code/upset_plot.R real/CRC/outcome/sig_family.csv CRC
```

This input table must be supplied separately; it is not automatically generated by the analysis scripts. The command saves `plot/real/CRC_overlap.pdf`. UpSet plotting writes only PDF files.

## Reproducibility notes

- Keep the directory structure intact; no absolute machine-specific paths need to be edited.
- Configure simulation workers before starting a full run. The default settings are intended for substantial parallel computation.
- Keep the original result CSV filenames so the bar plot script can identify lineage, target, signal, scenario, and inference type.
- Record `sessionInfo()` with your analysis results to document R and dependency versions. This repository does not include a dependency lockfile.
- Full simulations and real-data analyses must be run to obtain their statistical results; plotting scripts do not rerun these analyses.
