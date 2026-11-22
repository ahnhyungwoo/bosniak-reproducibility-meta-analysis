# Bosniak classification reproducibility meta-analysis

This repository contains the analytic data and R code for:

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
|   |   |-- management_tier_interreader.R
|   |   |-- management_tier_intermodality.R
|   |   |-- variance_swap_sensitivity.R
|   |   |-- reader_reference_sensitivity.R
|   |   |-- R2_full_table_models.R
|   |   `-- results/                 # R2 analysis tables
|   |-- figures/r2/                   # final-size figure data and renderer
|   `-- data/
|       |-- data_correction_log.csv
|       `-- version_matched_pair_audit.csv
|-- DATA_DICTIONARY.md
`-- run_all.R
```

The master extraction table contains 466 comparison records. Records with a
non-empty `exclude` field are retained for transparency but are not included in
the analysis. The final analytic dataset contains 332 comparisons from 79
studies.

Of these, 295 belong to the four principal strata, 16 are reader-versus-
reference comparisons (`reader-reference`), and 21 are technical or reading-
condition comparisons. The primary inter-reader set includes only independent
reader comparisons (159 effects, 46 studies, 51 clusters). For Münch, 84
examinations were read twice. The binary tables retain 168 paired ratings,
not 168 independent examinations.

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
checks for the publication-type contrasts in both primary strata.

The full runner establishes the repository root as its working directory.
Use it for the complete R2 results, including the coefficient-specific
management-tier fits, reader-reference add-back, table subgroups, and
publication-bias statistics. It does not require Python. To render Figures
2, 3A (two parts), 3B and S2A after the R run:

```sh
python -m pip install matplotlib Pillow
python revision/figures/r2/render_figures.py
```

Use Arial for the archived typography. The renderer creates 300-dpi masters,
72-dpi previews, a figure inventory and an audit of plotted comparison IDs.
The older plotting scripts in `meta_analysis/` are retained for provenance;
they do not generate the revised layouts.

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
  full-model validation checks, focused publication-type checks in both primary
  strata, management-tier analyses, and variance-source sensitivity analyses.
  `intermodality_management_tier_rows.csv` provides the per-comparison 2x2
  management-tier reconstructions underlying the inter-modality management-tier
  and directional synthesis reported in the Results.

See `DATA_DICTIONARY.md` for variable definitions and analytic conventions.

The audited 13 version pairs and their source locations are in
`revision/data/version_matched_pair_audit.csv`. In the results directory,
`version_matched_sensitivity.csv` contains the paired-correlation, inference,
and alternate reader-pair fits; `version_matched_leave_one_out.csv` deletes
each of the 12 comparable pairs in turn. `reader_reference_effects_by_study.csv`
identifies the 16 reference comparisons, and `reader_reference_sensitivity.csv`
refits the inter-reader model with and without them. `publication_bias_summary.csv`
contains Egger, Begg and trim-and-fill results for both OPES strata.
`R2_full_table_additional_models.csv` provides the remaining subgroup and
sensitivity rows; `primary_I2_percent.csv` separates between- and within-cluster
heterogeneity.

Two files under `results/` are not generated by `run_all.R`:
`intermodality_management_tier_rows.csv` is the extracted input to the
direction-of-crossing summaries; `version_all_effects_restricted.csv` is
a retained R1 ancillary result, not the R2 audited paired-version comparison.
Use `version_matched_summary.csv` and `version_matched_sensitivity.csv`
for that comparison.

## Analysis overview

The four principal agreement domains are inter-reader, inter-modality,
intra-reader, and inter-version agreement. Primary sufficiently populated
strata use multilevel random-effects models with cluster-robust variance
estimation and small-sample correction. Sparse strata use a one-per-study
estimate hierarchy. The primary three-level prediction interval targets the
true effect for a new individual comparison in a new cluster: it includes both
heterogeneity components. A separate interval for the latent cluster mean
includes only the between-cluster component.
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
