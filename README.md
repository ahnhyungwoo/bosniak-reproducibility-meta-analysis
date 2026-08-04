# Bosniak classification reproducibility meta-analysis

This repository contains the aggregate study-level data and R code for:

> Reproducibility of the Bosniak Classification for Malignancy Risk
> Stratification of Cystic Renal Masses: A Systematic Review and Meta-Analysis

## Repository contents

```text
.
|-- meta_analysis/
|   |-- comparisons.csv
|   |-- studies.csv
|   |-- qarel_for_analysis.csv
|   |-- 01_main_analysis.R
|   |-- 02_additional_sensitivity.R
|   |-- 03_publication_bias.R
|   |-- 04_forest_plots.R
|   `-- 04b_forest_compact.R
|-- revision/
|   |-- analysis/
|   |   |-- revision_analyses.R
|   |   `-- results/
|   `-- data/
|       `-- data_correction_log.csv
|-- DATA_DICTIONARY.md
`-- run_all.R
```

The master extraction table contains 466 comparison records. Records with a
non-empty `exclude` field are retained for transparency but are not included in
the analysis. The final analytic dataset contains 332 comparisons from 79
studies.

## Reproducing the analyses

The code was tested with R 4.4.2 and the following package versions:

- `metafor` 5.0-1
- `clubSandwich` 0.7.0
- `dplyr` 1.2.1
- `readr` 2.2.0

Install the required packages:

```r
install.packages(c("metafor", "clubSandwich", "dplyr", "readr"))
```

From the repository root, run:

```sh
Rscript run_all.R
```

To reproduce only the additional analyses performed during revision, run:

```sh
Rscript revision/analysis/revision_analyses.R
```

This command also reproduces the leave-one-study-out checks for all moderator
contrasts with nominal p<0.05 and the focused leave-one-abstract-study-out
checks for the inter-reader publication-type contrast.

Each script resolves its input files relative to its own location; no
author-specific working directory is required.

## Data files

- `comparisons.csv` contains one row per extracted agreement comparison,
  including the agreement statistic, sampling information, comparison type,
  dependency cluster, analytic inclusion status, and one-per-study estimate
  selection fields.
- `studies.csv` contains study-level bibliographic, design, population, imaging,
  and Bosniak-classification variables.
- `qarel_for_analysis.csv` contains the consensus assessments for seven
  QAREL-derived domains for the 79 studies retained in the analysis.
- `data_correction_log.csv` records corrections applied transparently during
  revision.
- `revision/analysis/results/` contains the tabular outputs from the revision
  analyses, including item-level QAREL counts and denominators, full
  leave-one-study-out iteration results, compact leave-one-study-out summaries,
  full-model validation checks, and the focused inter-reader publication-type
  checks.

See `DATA_DICTIONARY.md` for variable definitions and analytic conventions.

## Analysis overview

The four principal agreement domains are inter-reader, inter-modality,
intra-reader, and inter-version agreement. Primary sufficiently populated
strata use multilevel random-effects models with cluster-robust variance
estimation and small-sample correction. Sparse strata use a one-per-study
estimate hierarchy. Prediction intervals are reported where estimable.
Meta-regression is exploratory. Tests with Satterthwaite degrees of freedom
below 4, or a moderator level represented by only one study and one dependency
cluster, are treated as non-interpretable; those with degrees of freedom from 4
to below 10 are treated as fragile. Leave-one-study-out analyses are performed
for moderator contrasts meeting the protocol's nominal p<0.05 threshold, with
non-estimable deletion models and singleton moderator levels identified
explicitly in the output files.

## Citation

The manuscript is under revision. Until a final journal citation is available,
please cite this repository and the manuscript title shown above.

Repository: https://github.com/ahnhyungwoo/bosniak-reproducibility-meta-analysis
