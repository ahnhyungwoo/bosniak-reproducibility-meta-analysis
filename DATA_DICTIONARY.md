# Data dictionary

## `meta_analysis/comparisons.csv`

Each row is an extracted agreement comparison.

- `comparison_id`, `study_id`: stable comparison and study identifiers.
- `comparison_type`, `analytic_stratum`: extracted comparison type and the
  analysis stratum used in the meta-analysis.
- `modality_1`, `modality_2`, `std_modality_1`, `std_modality_2`: reported and
  standardized imaging modalities.
- `version_1`, `version_2`, `std_version_1`, `std_version_2`: reported and
  standardized Bosniak versions.
- `n_readers`, `n_cysts`, `n_categories`, `categories_used`,
  `collapsing_scheme`: sampling and category structure.
- `reader_1_specialty`, `reader_2_specialty`,
  `reader_1_experience_yrs`, `reader_2_experience_yrs`, `reader_structure`,
  `specialty_group`, `exp_mean`: reader characteristics and derived groupings.
- `blinding`, `blinding_reported`, `reading_session_type`: reported reading
  conditions.
- `kappa`, `gwet_ac1`, `gwet_ac2`, `icc`, `krippendorff_alpha`: extracted
  agreement statistics.
- `kappa_type`, `weight_scheme`, `icc_type`: statistic specifications.
- `kappa_se`, `ci_lower`, `ci_upper`, `observed_agreement`: uncertainty and
  observed agreement, when available.
- `has_raw_table`, `metric_note`: whether a cross-tabulation was available and
  relevant metric provenance.
- `exclude`: non-empty values identify records excluded from the analysis and
  give the exclusion reason. The value `not_inter_reader_agreement` identifies
  reader-specific diagnostic-agreement values that are not inter-reader
  reliability estimates.
- `variation_axis`, `sensitivity_flags`: comparison structure and sensitivity
  analysis flags.
- `dep_id`: dependency cluster used for robust variance estimation.
- `opes_rank`, `opes_include`, `opes_basis`, `opes_exclude_reason`: fields for
  selecting one representative estimate per study in sparse strata.
- `study_design_cat`, `region`: standardized study-level moderator fields.

Blank cells represent information that was not reported or was not applicable.

## `meta_analysis/studies.csv`

Each row is a study-level record linked to `comparisons.csv` by `study_id`.
Fields include bibliographic metadata, publication type, country, study design,
enrollment period, setting, population selection, sample size, age and sex,
lesion size, Bosniak-category counts, malignancy rate, reference standard,
Bosniak version, and imaging acquisition characteristics.

## `meta_analysis/qarel_for_analysis.csv`

The file contains one row for each of the 79 studies retained in the analysis
and is linked to the other data files by `study_id`.

- `Q1`-`Q7`: consensus ratings for the seven QAREL-derived domains used in the
  review.
- `qarel_yes_count`, `qarel_total`: number of affirmative ratings and the
  number of assessed domains.
- `overall_rob`: derived overall risk-of-bias category.

`revision/analysis/results/qarel_item_counts.csv` reports Yes, No, Unclear,
Not Applicable, and missing counts for each domain, together with the
all-study and applicable-study denominators and percentages.

## Revision analysis outputs

- `loso_nominal_moderator_iterations.csv` contains one row for every planned
  leave-one-study-out deletion for the eight moderator contrasts with nominal
  p<0.05. It records model estimability, the omitted study, the coefficient,
  standard error, Satterthwaite degrees of freedom, raw p value, coefficient
  direction, and warnings or failure reason.
- `loso_nominal_moderator_summary.csv` summarizes the successful and
  non-estimable deletion models, direction consistency, nominal p<0.05 counts,
  and coefficient and p-value ranges for each of the eight contrasts.
- `loso_nominal_moderator_validation.csv` verifies that every recomputed full
  model matches its row in `meta_regression_revision.csv` before the deletion
  analyses are run.
- `inter_reader_publication_type_leave_one_abstract.csv` contains the full
  inter-reader publication-type model and the two focused checks omitting each
  abstract study in turn, with singleton abstract levels identified as
  non-interpretable.

## Analytic conventions

- Records are analytically active when `exclude` is empty.
- Agreement coefficients are clamped to the interval -0.995 to 0.995 before
  Fisher transformation.
- Sampling variances are taken from reported confidence intervals or standard
  errors when available; otherwise they are approximated from the lesion count.
- Linear, quadratic, custom, Cicchetti, and weighted-but-unspecified kappa
  estimates are grouped as weighted for the relevant sensitivity analyses.
- Primary sufficiently populated strata use multilevel random-effects models
  with cluster-robust variance estimation and CR2 small-sample correction.
- Sparse strata use the `opes_*` hierarchy to select one representative
  estimate per study. The analysis stops if any eligible study-stratum has
  zero or more than one selected estimate.
