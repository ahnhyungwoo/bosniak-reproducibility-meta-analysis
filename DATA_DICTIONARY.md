# Data dictionary

## `meta_analysis/comparisons.csv`

Each row is an extracted agreement comparison.

- `comparison_id`, `study_id`: stable comparison and study identifiers.
- `comparison_type`, `analytic_stratum`: extracted comparison type and the
  analysis stratum used in the meta-analysis.
  `reader-reference` distinguishes retained comparisons against consensus or
  an expert standard from independent inter-reader comparisons. It is neither
  a fifth primary stratum nor an exclusion flag. Sixteen such comparisons
  remain in the denominator of 332, but not the 203 direct inter-reader records.
- `modality_1`, `modality_2`, `std_modality_1`, `std_modality_2`: reported and
  standardized imaging modalities.
- `version_1`, `version_2`, `std_version_1`, `std_version_2`: reported and
  standardized Bosniak versions.
- `n_readers`, `n_cysts`, `n_categories`, `categories_used`,
  `collapsing_scheme`: sampling and category structure.
  For Münch's three direct inter-reader comparisons and six comparisons
  against consensus, `n_cysts=84` denotes independently evaluated CEUS
  examinations; 168 is the number of presentations to each reader across two
  rounds. Reported CIs take precedence over the sample-size proxy.
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

### `intermodality_management_tier_rows.csv`

This file contains the 26 CT-MRI and CT-CEUS management-tier comparison rows
underlying the inter-modality management-tier and directional synthesis reported
in the Results. Classifications were collapsed at I-IIF versus III-IV or the
source-defined equivalent, with CT as modality 1 and MRI or CEUS as modality 2.

- `tier_row_id`: stable public row identifier (`IMT01`-`IMT26`).
- `modality_pair`, `study`, `version_reader`, `n`: comparison context and lesion
  count.
- `comparison_id`: link to `meta_analysis/comparisons.csv` where a unique
  effect-level record was available; blank for nine reconstructed CT-CEUS rows.
- `lower_lower`: CT lower tier and modality 2 lower tier.
- `other_higher`: CT lower tier and modality 2 higher tier (an upward crossing).
- `other_lower`: CT higher tier and modality 2 lower tier (a downward crossing).
- `higher_higher`: CT higher tier and modality 2 higher tier.
- `observed_agreement`, `kappa`, `ac1`: per-comparison agreement estimates.
- `net_other_higher`: `other_higher - other_lower`; positive, zero, and negative
  values define upward, balanced, and downward net tier shifts, respectively.
- `reconstruction_status`: `full` when all four 2x2 cells were reconstructable
  and `direction-only` when only directional crossings and observed agreement
  were available.
- `source`: source-table or report location used for reconstruction.

For Rocca (2016), `lower_lower` and `higher_higher` are blank because the
directional crossings and observed agreement did not identify how concordant
lesions were divided between the lower and higher tiers; marginal distributions,
κ, and AC1 therefore could not be reconstructed.

### Management-tier analysis outputs

- `interreader_management_tier_cells.csv` contains the reconstructed lower- and
  higher-tier 2x2 cells and derived agreement measures for every assessed
  inter-reader comparison.
- `interreader_management_tier_summary.csv` contains the RVE+CR2 pooled κ and
  AC1 analyses at the management-tier boundary.
- `interreader_management_tier_descriptive.csv` reports comparison counts,
  study counts, and observed-agreement ranges and medians for the reconstructed
  inter-reader sets.
- `interreader_management_tier_accounting.csv` records the reconstruction and
  inclusion accounting, including excluded or non-estimable tables.
- `intermodality_management_tier_summary.csv` summarizes the full CT-MRI and
  CT-CEUS reconstructions and the direction-only Rocca row from
  `intermodality_management_tier_rows.csv`.

### Additional sensitivity outputs

- `variance_swap_sensitivity.csv` reports the post hoc sampling-variance
  substitution sensitivity analysis added during revision.
- `version_all_effects_restricted.csv` reports the restricted all-effects
  version sensitivity analyses.

- `loso_nominal_moderator_iterations.csv` contains one row for every planned
  leave-one-study-out deletion for the six moderator contrasts with nominal
  p<0.05. It records model estimability, the omitted study, the coefficient,
  standard error, Satterthwaite degrees of freedom, raw p value, coefficient
  direction, and warnings or failure reason.
- `loso_nominal_moderator_summary.csv` summarizes the successful and
  non-estimable deletion models, direction consistency, nominal p<0.05 counts,
  and coefficient and p-value ranges for each of the six contrasts.
- `loso_nominal_moderator_validation.csv` verifies that every recomputed full
  model matches its row in `meta_regression_revision.csv` before the deletion
  analyses are run.
- `publication_type_leave_one_abstract.csv` contains the full publication-type
  models for the inter-reader and inter-modality strata and the focused checks
  omitting each abstract study in turn, with singleton abstract levels
  identified as non-interpretable.
- `inter_reader_publication_type_leave_one_abstract.csv` is retained as an
  inter-reader-only extract for compatibility.
- `possible_cohort_overlap_sensitivity.csv` contains the inter-reader
  single-cohort dependency stress test for Bai (2020) and Kang (2022), the
  corresponding full OPES estimate, and OPES estimates excluding each study in
  turn.

## Analytic conventions

- `revision/data/version_matched_pair_audit.csv`: source file/location, matched
  comparison IDs, reader/lesion/scale/weighting checks and effective counts
  for each version. A common cohort is not always the same evaluable set.
  `version_matched_pairs.csv` carries the actual effect and variance inputs.
  `version_matched_sensitivity.csv`: 12 comparable pairs, all 13 pairs, and
  alternative Yan pairs; within-pair correlations and normal/Knapp–Hartung
  inference. The main estimand is the directly pooled within-study difference
  in κ (v2019 minus pre-2019), not a difference between pooled levels.
  `version_matched_leave_one_out.csv`: excluded study and its κ difference,
  remaining pair count, pooled κ difference, CI and p (REML, KH, r=0).
- `reader_reference_effects_by_study.csv`: study, IDs and descriptive κ
  range/median; `reader_reference_sensitivity.csv` distinguishes the primary
  direct-reader model from the 16-effect add-back. The original estimates
  can be retrieved by comparison ID.
- In `interreader_management_tier_cells.csv`, `n` is the observed 2×2
  cell total for an included table; `independent_lesions` separately records
  84 for MTIR031–038 (blank elsewhere means not separately audited).
  For these eight tables only, the proxy variance `vi` and the two Fisher-z
  scale jackknife variances `vi_kappa_jk` and `vi_ac1_jk` are doubled as a
  conservative approximation for repeated readings, not a lesion-level
  jackknife. Coefficients and cell counts are unchanged.
- In three-level output, `sigma2_between` and `sigma2_within` are variance
  components on the model scale. `pi_lower/upper` predict a new comparison
  in a new cluster and include both; `pi_cluster_mean_lower/upper` predict a
  latent cluster mean and include only the between-cluster variance. Both
  include estimation uncertainty. Single-level OPES intervals have no
  separate comparison-level component.
- `publication_bias_summary.csv`: one row per OPES stratum, study count,
  Egger regression t/df/p, Kendall tau and p, imputed study count and side,
  pooled and trim-and-fill κ and CIs. These assess small-study effects, not
  a diagnosis of publication bias.

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
