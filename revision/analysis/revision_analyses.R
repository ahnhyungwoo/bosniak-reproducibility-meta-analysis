## Bosniak reproducibility major-revision analyses
##
## This script leaves the submitted datasets unchanged and writes the
## additional analysis outputs under revision/analysis/results.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run this file with Rscript.")
}
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
project_dir <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/")
local_lib <- file.path(project_dir, "revision", "R_lib")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(metafor)
  library(clubSandwich)
  library(dplyr)
  library(readr)
})

data_dir <- file.path(project_dir, "meta_analysis")
out_dir <- file.path(project_dir, "revision", "analysis", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

comp <- read_csv(file.path(data_dir, "comparisons.csv"), show_col_types = FALSE)
stud <- read_csv(file.path(data_dir, "studies.csv"), show_col_types = FALSE)
qarel <- read_csv(file.path(data_dir, "qarel_for_analysis.csv"), show_col_types = FALSE)

active_study_ids <- comp %>%
  filter(is.na(exclude) | trimws(exclude) == "") %>%
  distinct(study_id)

qarel_analytic <- qarel %>%
  semi_join(active_study_ids, by = "study_id")

if (n_distinct(qarel_analytic$study_id) != nrow(active_study_ids)) {
  stop("QAREL records do not match the retained-study set.")
}

qarel_item_labels <- c(
  Q1 = "Representative spectrum",
  Q2 = "Adequate rater description",
  Q3 = "Independent or blinded interpretation",
  Q4 = "Blinding to clinical information",
  Q5 = "Appropriate time interval or washout",
  Q6 = "Classification criteria specified",
  Q7 = "Appropriate statistical analysis"
)

qarel_item_counts <- bind_rows(lapply(names(qarel_item_labels), function(item) {
  rating <- qarel_analytic[[item]]
  not_applicable <- rating %in% c("N/A", "Not Applicable")
  applicable <- !is.na(rating) & trimws(rating) != "" & !not_applicable
  yes_count <- sum(rating == "Yes", na.rm = TRUE)

  tibble(
    item = item,
    domain = unname(qarel_item_labels[[item]]),
    studies_total = nrow(qarel_analytic),
    yes = yes_count,
    no = sum(rating == "No", na.rm = TRUE),
    unclear = sum(rating == "Unclear", na.rm = TRUE),
    not_applicable = sum(not_applicable, na.rm = TRUE),
    missing = sum(is.na(rating) | trimws(rating) == ""),
    applicable_denominator = sum(applicable),
    yes_pct_all_studies = round(100 * yes_count / nrow(qarel_analytic), 1),
    yes_pct_applicable = round(100 * yes_count / sum(applicable), 1)
  )
}))

write_csv(
  qarel_item_counts,
  file.path(out_dir, "qarel_item_counts.csv")
)

# Source verification during revision showed that Chang 2015 is a poster
# abstract in the 2015 JASN conference supplement, not a full journal article.
# The submitted source file is left untouched; the correction is applied
# transparently here and documented in revision/data/data_correction_log.csv.
stud <- stud %>%
  mutate(
    publication_type = if_else(
      study_id == "Chang_2015",
      "conference_abstract",
      publication_type
    )
  )

pediatric_ids <- c("Peng_2010", "Karmazyn_2015", "Frumer_2021", "Peard_2022")

dat <- comp %>%
  filter(is.na(exclude) | trimws(exclude) == "") %>%
  left_join(
    stud %>%
      select(
        study_id, year, publication_type, study_design, bosniak_version,
        n_patients_study = n_patients, n_cysts_study = n_cysts
      ),
    by = "study_id"
  ) %>%
  left_join(
    qarel %>% select(study_id, Q1, overall_rob, qarel_yes_count),
    by = "study_id"
  ) %>%
  mutate(
    kappa_num = suppressWarnings(as.numeric(kappa)),
    ci_lo = suppressWarnings(as.numeric(ci_lower)),
    ci_hi = suppressWarnings(as.numeric(ci_upper)),
    se_raw = suppressWarnings(as.numeric(kappa_se)),
    nc = suppressWarnings(as.numeric(n_cysts)),
    n_readers_num = suppressWarnings(as.numeric(n_readers)),
    n_categories_num = suppressWarnings(as.numeric(n_categories)),
    exp_mean_num = suppressWarnings(as.numeric(exp_mean)),
    year_num = suppressWarnings(as.numeric(year)),
    kappa_clamped = pmin(pmax(kappa_num, -0.995), 0.995),
    yi_z = atanh(kappa_clamped),
    variance_source = case_when(
      !is.na(ci_lo) & !is.na(ci_hi) ~ "reported_CI",
      !is.na(se_raw) ~ "reported_SE",
      !is.na(nc) & nc > 3 ~ "n_approximation",
      TRUE ~ "unavailable"
    ),
    vi_z_ci = if_else(
      !is.na(ci_lo) & !is.na(ci_hi),
      ((atanh(pmin(pmax(ci_hi, -0.995), 0.995)) -
          atanh(pmin(pmax(ci_lo, -0.995), 0.995))) /
         (2 * qnorm(0.975)))^2,
      NA_real_
    ),
    vi_z_se = if_else(
      !is.na(se_raw),
      (se_raw / (1 - kappa_clamped^2))^2,
      NA_real_
    ),
    vi_z_n = if_else(!is.na(nc) & nc > 3, 1 / (nc - 3), NA_real_),
    vi_z = coalesce(vi_z_ci, vi_z_se, vi_z_n),
    vi_k_ci = if_else(
      !is.na(ci_lo) & !is.na(ci_hi),
      ((ci_hi - ci_lo) / (2 * qnorm(0.975)))^2,
      NA_real_
    ),
    vi_k_se = if_else(!is.na(se_raw), se_raw^2, NA_real_),
    vi_k_n = if_else(
      !is.na(nc) & nc > 3,
      (1 - kappa_clamped^2)^2 / (nc - 3),
      NA_real_
    ),
    vi_kappa = coalesce(vi_k_ci, vi_k_se, vi_k_n),
    wt_group = case_when(
      weight_scheme %in%
        c("linear", "quadratic", "custom", "Cicchetti", "weighted_unspecified") ~
        "weighted",
      weight_scheme == "unweighted" ~ "unweighted",
      TRUE ~ "unspecified"
    ),
    ncat_group = case_when(
      n_categories_num <= 3 ~ "2-3",
      n_categories_num <= 5 ~ "4-5",
      n_categories_num > 5 ~ "6+",
      TRUE ~ NA_character_
    ),
    mod_group = case_when(
      std_modality_1 == "CT" ~ "CT",
      std_modality_1 == "MRI" ~ "MRI",
      std_modality_1 %in% c("CEUS", "US", "SMI") ~ "CEUS_US",
      TRUE ~ "other"
    ),
    mod_pair = case_when(
      (std_modality_1 %in% c("CT", "DECT") & std_modality_2 == "MRI") |
        (std_modality_1 == "MRI" & std_modality_2 %in% c("CT", "DECT")) ~
        "CT_MRI",
      (std_modality_1 %in% c("CT", "CT_MRI") &
         std_modality_2 %in% c("CEUS", "US")) |
        (std_modality_1 %in% c("CEUS", "US") &
           std_modality_2 %in% c("CT", "CT_MRI")) ~
        "CT_CEUS_US",
      (std_modality_1 == "MRI" &
         std_modality_2 %in% c("CEUS", "US")) |
        (std_modality_1 %in% c("CEUS", "US") &
           std_modality_2 == "MRI") ~
        "MRI_CEUS_US",
      TRUE ~ "other"
    ),
    pub_group = case_when(
      publication_type == "conference_abstract" ~ "abstract",
      publication_type == "journal_article" ~ "fulltext",
      TRUE ~ NA_character_
    ),
    log_nc = if_else(!is.na(nc) & nc > 0, log(nc), NA_real_),
    spec_group = case_when(
      specialty_group == "subspecialty" ~ "subspecialty",
      specialty_group %in%
        c("general", "mixed", "residents", "non_radiology") ~
        "non_subspecialty",
      TRUE ~ NA_character_
    ),
    blind_group = case_when(
      blinding_reported == "blinded" ~ "reported",
      blinding_reported %in% c("unblinded", "unclear") ~ "not_reported",
      TRUE ~ NA_character_
    ),
    structure_group = case_when(
      reader_structure == "pooled" ~ "pooled",
      reader_structure %in%
        c("single_pair", "pairwise", "single", "subgroup") ~
        "other",
      TRUE ~ NA_character_
    ),
    design_group = case_when(
      study_design_cat == "prospective" ~ "prospective",
      study_design_cat == "retrospective" ~ "retrospective",
      TRUE ~ NA_character_
    ),
    region_group = if_else(
      region %in% c("North_America", "Asia", "Europe"),
      region,
      NA_character_
    ),
    representative_spectrum = if_else(Q1 == "Yes", "Yes", "No_or_unclear"),
    standard_v1 = std_version_1 %in% c("original", "v2019"),
    standard_v2 = std_version_2 %in% c("original", "v2019"),
    standard_adult = !(study_id %in% pediatric_ids) &
      case_when(
        analytic_stratum == "inter-modality" ~ standard_v1 & standard_v2,
        TRUE ~ standard_v1
      )
  ) %>%
  filter(!is.na(kappa_num), !is.na(vi_z), vi_z > 0)

ir <- dat %>% filter(analytic_stratum == "inter-reader")
im <- dat %>% filter(analytic_stratum == "inter-modality")
intra <- dat %>% filter(analytic_stratum == "intra-reader")
iv <- dat %>% filter(analytic_stratum == "inter-version")

validate_opes_selection <- function(data) {
  audit <- data %>%
    filter(analytic_stratum %in% c(
      "inter-reader", "inter-modality", "intra-reader", "inter-version"
    )) %>%
    group_by(analytic_stratum, study_id) %>%
    summarise(
      selected_effects = sum(opes_include == 1, na.rm = TRUE),
      selected_clusters = n_distinct(
        dep_id[coalesce(opes_include == 1, FALSE)], na.rm = TRUE
      ),
      .groups = "drop"
    )

  invalid <- audit %>%
    filter(selected_effects != 1 | selected_clusters != 1)
  if (nrow(invalid) > 0) {
    details <- paste(
      sprintf(
        "%s/%s: effects=%d, clusters=%d",
        invalid$analytic_stratum,
        invalid$study_id,
        invalid$selected_effects,
        invalid$selected_clusters
      ),
      collapse = "; "
    )
    stop("Invalid one-per-study estimate selection: ", details)
  }
}

validate_opes_selection(dat)

region_levels <- c("North America", "Asia", "Europe", "Other/unspecified region")
ir_region_data <- ir %>%
  mutate(
    region_display = case_when(
      region == "North_America" ~ "North America",
      region == "Asia" ~ "Asia",
      region == "Europe" ~ "Europe",
      TRUE ~ "Other/unspecified region"
    ),
    version_display = if_else(
      !is.na(std_version_1) & std_version_1 == "v2019",
      "2019 update",
      "Non-2019"
    ),
    modality_display = case_when(
      std_modality_1 == "CT" ~ "CT",
      std_modality_1 == "MRI" ~ "MRI",
      std_modality_1 %in% c("CEUS", "US", "SMI") ~ "CEUS/US",
      std_modality_1 == "CT_MRI" ~ "CT/MRI",
      TRUE ~ "Mixed/not reported"
    ),
    publication_display = if_else(
      publication_type == "conference_abstract",
      "Conference abstract",
      "Journal article"
    )
  )

complete_region_counts <- function(data, section, variable, level_order) {
  observed <- data %>%
    count(region_display, level = .data[[variable]], name = "count")
  grid <- as_tibble(expand.grid(
    region_display = region_levels,
    level = level_order,
    stringsAsFactors = FALSE
  ))
  totals <- data %>% count(region_display, name = "region_effects")

  grid %>%
    left_join(observed, by = c("region_display", "level")) %>%
    left_join(totals, by = "region_display") %>%
    mutate(
      section = section,
      count = coalesce(count, 0L),
      percent = round(100 * count / region_effects),
      region_order = match(region_display, region_levels),
      level_order = match(level, level_order)
    ) %>%
    arrange(region_order, level_order) %>%
    select(section, level, region_display, region_effects, count, percent)
}

inter_reader_region_distribution <- bind_rows(
  complete_region_counts(
    ir_region_data,
    "Bosniak version",
    "version_display",
    c("Non-2019", "2019 update")
  ),
  complete_region_counts(
    ir_region_data,
    "Imaging modality",
    "modality_display",
    c("CT", "MRI", "CEUS/US", "CT/MRI", "Mixed/not reported")
  ),
  complete_region_counts(
    ir_region_data,
    "Publication type",
    "publication_display",
    c("Journal article", "Conference abstract")
  )
)

write_csv(
  inter_reader_region_distribution,
  file.path(out_dir, "inter_reader_region_distribution.csv")
)

safe_tcrit <- function(df) {
  ifelse(is.finite(df) & df > 0, qt(0.975, df), qnorm(0.975))
}

run_rve <- function(
    data,
    label,
    outcome = c("z", "kappa"),
    cluster_mode = c("dependency", "study"),
    analysis_family = "sensitivity") {
  outcome <- match.arg(outcome)
  cluster_mode <- match.arg(cluster_mode)
  if (nrow(data) < 3) return(NULL)

  d <- data
  d$analysis_yi <- if (outcome == "z") d$yi_z else d$kappa_num
  d$analysis_vi <- if (outcome == "z") d$vi_z else d$vi_kappa
  d <- d %>% filter(!is.na(analysis_yi), !is.na(analysis_vi), analysis_vi > 0)
  if (nrow(d) < 3) return(NULL)

  if (cluster_mode == "dependency") {
    cluster_var <- d$dep_id
    random_spec <- list(~ 1 | dep_id, ~ 1 | comparison_id)
  } else {
    cluster_var <- d$study_id
    random_spec <- list(~ 1 | study_id, ~ 1 | comparison_id)
  }
  if (length(unique(cluster_var)) < 3) return(NULL)

  model <- rma.mv(
    analysis_yi,
    analysis_vi,
    random = random_spec,
    data = d,
    method = "REML"
  )
  robust <- coef_test(model, vcov = "CR2", cluster = cluster_var)
  estimate_working <- as.numeric(model$b[1])
  se <- robust$SE[1]
  df <- robust$df_Satt[1]
  crit <- safe_tcrit(df)
  ci_working <- estimate_working + c(-1, 1) * crit * se
  pi_working <- estimate_working + c(-1, 1) * crit *
    sqrt(se^2 + model$sigma2[1])

  if (outcome == "z") {
    estimate <- tanh(estimate_working)
    ci <- tanh(ci_working)
    pi <- tanh(pi_working)
  } else {
    estimate <- estimate_working
    ci <- pmin(pmax(ci_working, -1), 1)
    pi <- pmin(pmax(pi_working, -1), 1)
  }

  tibble(
    analysis_family = analysis_family,
    analysis = label,
    outcome_scale = outcome,
    cluster_mode = cluster_mode,
    effects = nrow(d),
    clusters = length(unique(cluster_var)),
    studies = n_distinct(d$study_id),
    estimate = estimate,
    ci_lower = ci[1],
    ci_upper = ci[2],
    pi_lower = pi[1],
    pi_upper = pi[2],
    robust_se_working_scale = se,
    satterthwaite_df = df,
    p_value = robust$p_Satt[1]
  )
}

prepare_opes <- function(data) {
  data %>%
    filter(opes_include == 1) %>%
    mutate(
      opes_value_text = if_else(
        grepl("^(median_synthetic|reported_overall_restored):", opes_basis),
        sub("^[^:]+:", "", opes_basis),
        NA_character_
      ),
      opes_kappa = coalesce(suppressWarnings(as.numeric(opes_value_text)), kappa_num),
      opes_kappa = pmin(pmax(opes_kappa, -0.995), 0.995),
      opes_yi = atanh(opes_kappa),
      opes_vi = if_else(
        !is.na(opes_value_text),
        if_else(!is.na(nc) & nc > 3, 1 / (nc - 3), vi_z),
        vi_z
      )
    ) %>%
    filter(!is.na(opes_yi), !is.na(opes_vi), opes_vi > 0)
}

run_opes <- function(data, label, analysis_family = "primary") {
  d <- prepare_opes(data)
  if (nrow(d) < 2) return(NULL)
  model <- rma(opes_yi, opes_vi, data = d, method = "REML")
  pred <- predict(model)
  tibble(
    analysis_family = analysis_family,
    analysis = label,
    outcome_scale = "z",
    cluster_mode = "one_per_study",
    effects = nrow(d),
    clusters = n_distinct(d$study_id),
    studies = n_distinct(d$study_id),
    estimate = tanh(as.numeric(model$b[1])),
    ci_lower = tanh(model$ci.lb),
    ci_upper = tanh(model$ci.ub),
    pi_lower = tanh(pred$pi.lb),
    pi_upper = tanh(pred$pi.ub),
    robust_se_working_scale = model$se,
    satterthwaite_df = NA_real_,
    p_value = model$pval
  )
}

run_rho <- function(data, label, rho) {
  d <- data %>% arrange(dep_id, comparison_id)
  if (nrow(d) < 3 || n_distinct(d$dep_id) < 3) return(NULL)
  working_v <- impute_covariance_matrix(
    vi = d$vi_z,
    cluster = d$dep_id,
    r = rho,
    smooth_vi = FALSE,
    return_list = FALSE
  )
  model <- rma.mv(
    yi_z,
    V = working_v,
    random = list(~ 1 | dep_id, ~ 1 | comparison_id),
    data = d,
    method = "REML"
  )
  robust <- coef_test(model, vcov = "CR2", cluster = d$dep_id)
  crit <- safe_tcrit(robust$df_Satt[1])
  ci_z <- as.numeric(model$b[1]) + c(-1, 1) * crit * robust$SE[1]
  tibble(
    analysis_family = "rho_sensitivity",
    analysis = paste0(label, " (rho=", rho, ")"),
    outcome_scale = "z",
    cluster_mode = "dependency_correlated_working_V",
    effects = nrow(d),
    clusters = n_distinct(d$dep_id),
    studies = n_distinct(d$study_id),
    estimate = tanh(as.numeric(model$b[1])),
    ci_lower = tanh(ci_z[1]),
    ci_upper = tanh(ci_z[2]),
    pi_lower = NA_real_,
    pi_upper = NA_real_,
    robust_se_working_scale = robust$SE[1],
    satterthwaite_df = robust$df_Satt[1],
    p_value = robust$p_Satt[1]
  )
}

summary_parts <- list(
  run_rve(ir, "Inter-reader primary RVE+CR2", analysis_family = "primary"),
  run_rve(im, "Inter-modality primary RVE+CR2", analysis_family = "primary"),
  run_opes(intra, "Intra-reader primary OPES", analysis_family = "primary"),
  run_opes(iv, "Inter-version primary OPES", analysis_family = "primary"),
  run_rve(
    ir %>% filter(standard_adult),
    "Inter-reader standard-scheme adult-only",
    analysis_family = "standard_adult"
  ),
  run_rve(
    im %>% filter(standard_adult),
    "Inter-modality standard-scheme adult-only",
    analysis_family = "standard_adult"
  ),
  run_opes(
    im %>% filter(standard_adult, mod_pair == "CT_CEUS_US"),
    "CT-US-based standard-scheme adult-only OPES",
    analysis_family = "standard_adult"
  ),
  run_rve(
    ir %>% filter(variance_source %in% c("reported_CI", "reported_SE")),
    "Inter-reader directly reported variance only",
    analysis_family = "variance_source"
  ),
  run_rve(
    im %>% filter(variance_source %in% c("reported_CI", "reported_SE")),
    "Inter-modality directly reported variance only",
    analysis_family = "variance_source"
  ),
  run_rve(
    ir,
    "Inter-reader untransformed kappa",
    outcome = "kappa",
    analysis_family = "transformation"
  ),
  run_rve(
    im,
    "Inter-modality untransformed kappa",
    outcome = "kappa",
    analysis_family = "transformation"
  ),
  run_rve(
    ir,
    "Inter-reader study-level clustering",
    cluster_mode = "study",
    analysis_family = "cluster"
  ),
  run_rve(
    im,
    "Inter-modality study-level clustering",
    cluster_mode = "study",
    analysis_family = "cluster"
  ),
  run_rve(
    ir %>% filter(nc >= 30),
    "Inter-reader n>=30 lesions",
    analysis_family = "sample_size"
  ),
  run_rve(
    im %>% filter(nc >= 30),
    "Inter-modality n>=30 lesions",
    analysis_family = "sample_size"
  ),
  run_rve(
    ir %>% filter(representative_spectrum == "Yes"),
    "Inter-reader representative-spectrum studies",
    analysis_family = "representative_spectrum"
  ),
  run_rve(
    im %>% filter(representative_spectrum == "Yes"),
    "Inter-modality representative-spectrum studies",
    analysis_family = "representative_spectrum"
  )
)

for (rho in c(0, 0.3, 0.5, 0.8)) {
  summary_parts <- append(
    summary_parts,
    list(
      run_rho(ir, "Inter-reader", rho),
      run_rho(im, "Inter-modality", rho)
    )
  )
}

analysis_summary <- bind_rows(summary_parts)
write_csv(analysis_summary, file.path(out_dir, "revision_analysis_summary.csv"))

variance_counts <- dat %>%
  filter(analytic_stratum %in% c(
    "inter-reader", "inter-modality", "intra-reader", "inter-version"
  )) %>%
  group_by(analytic_stratum, variance_source) %>%
  summarise(
    effects = n(),
    studies = n_distinct(study_id),
    clusters = n_distinct(dep_id),
    .groups = "drop"
  )
write_csv(variance_counts, file.path(out_dir, "variance_source_counts.csv"))

weight_counts <- dat %>%
  filter(analytic_stratum %in% c("inter-reader", "inter-modality")) %>%
  mutate(weight_scheme_report = coalesce(weight_scheme, "unspecified")) %>%
  group_by(analytic_stratum, weight_scheme_report) %>%
  summarise(
    effects = n(),
    studies = n_distinct(study_id),
    clusters = n_distinct(dep_id),
    .groups = "drop"
  )
write_csv(weight_counts, file.path(out_dir, "weight_scheme_counts.csv"))

active_comp <- comp %>%
  filter(is.na(exclude) | trimws(exclude) == "")
primary_strata <- c(
  "inter-reader", "inter-modality", "intra-reader", "inter-version"
)
included_study_ids <- sort(unique(active_comp$study_id))

comparison_accounting <- tibble(
  item = c(
    "Included studies",
    "All retained comparisons",
    "Primary-stratum comparisons",
    "Inter-reader comparisons",
    "Inter-modality comparisons",
    "Intra-reader comparisons",
    "Inter-version comparisons",
    "Technical/read-condition comparisons",
    "Studies reporting patient count",
    "Patients across reporting studies",
    "Studies reporting cyst count",
    "Cysts across reporting studies"
  ),
  value = c(
    length(included_study_ids),
    nrow(active_comp),
    sum(active_comp$analytic_stratum %in% primary_strata),
    sum(active_comp$analytic_stratum == "inter-reader"),
    sum(active_comp$analytic_stratum == "inter-modality"),
    sum(active_comp$analytic_stratum == "intra-reader"),
    sum(active_comp$analytic_stratum == "inter-version"),
    sum(active_comp$analytic_stratum == "technical-comparison"),
    sum(
      stud$study_id %in% included_study_ids &
        !is.na(suppressWarnings(as.numeric(stud$n_patients)))
    ),
    sum(
      suppressWarnings(as.numeric(
        stud$n_patients[stud$study_id %in% included_study_ids]
      )),
      na.rm = TRUE
    ),
    sum(
      stud$study_id %in% included_study_ids &
        !is.na(suppressWarnings(as.numeric(stud$n_cysts)))
    ),
    sum(
      suppressWarnings(as.numeric(
        stud$n_cysts[stud$study_id %in% included_study_ids]
      )),
      na.rm = TRUE
    )
  )
)
write_csv(
  comparison_accounting,
  file.path(out_dir, "comparison_accounting.csv")
)

reported_by_stratum <- active_comp %>%
  filter(analytic_stratum %in% primary_strata) %>%
  mutate(
    reported_metric = case_when(
      !is.na(suppressWarnings(as.numeric(kappa))) ~ "kappa",
      !is.na(observed_agreement) & trimws(observed_agreement) != "" ~
        "percentage agreement",
      !is.na(gwet_ac1) & trimws(gwet_ac1) != "" ~ "Gwet AC1",
      !is.na(gwet_ac2) & trimws(gwet_ac2) != "" ~ "Gwet AC2",
      !is.na(icc) & trimws(icc) != "" ~ "ICC",
      !is.na(krippendorff_alpha) & trimws(krippendorff_alpha) != "" ~
        "Krippendorff alpha",
      TRUE ~ "other/non-kappa"
    )
  ) %>%
  group_by(study_id, analytic_stratum) %>%
  summarise(
    retained_comparisons = n(),
    retained_metric_types = paste(sort(unique(reported_metric)), collapse = "; "),
    .groups = "drop"
  )

kappa_by_stratum <- dat %>%
  filter(analytic_stratum %in% primary_strata) %>%
  group_by(study_id, analytic_stratum) %>%
  summarise(
    eligible_kappa_effects = n(),
    selected_opes_effects = sum(opes_include == 1, na.rm = TRUE),
    .groups = "drop"
  )

study_stratum_accounting <- reported_by_stratum %>%
  left_join(
    stud %>%
      select(study_id, first_author, year, publication_type),
    by = "study_id"
  ) %>%
  left_join(
    kappa_by_stratum,
    by = c("study_id", "analytic_stratum")
  ) %>%
  mutate(
    eligible_kappa_effects = coalesce(eligible_kappa_effects, 0L),
    selected_opes_effects = coalesce(selected_opes_effects, 0L),
    primary_model_status = case_when(
      analytic_stratum %in% c("inter-reader", "inter-modality") &
        eligible_kappa_effects > 0 ~
        "Included in primary all-eligible-kappa RVE+CR2 model",
      analytic_stratum %in% c("inter-reader", "inter-modality") ~
        "Excluded from kappa pooling: only non-kappa metric(s) reported",
      analytic_stratum %in% c("intra-reader", "inter-version") &
        selected_opes_effects > 0 ~
        "Included in primary OPES model (one selected effect per study)",
      analytic_stratum %in% c("intra-reader", "inter-version") &
        eligible_kappa_effects > 0 ~
        paste(
          "Eligible dependent kappa effect(s) retained in all-effects",
          "sensitivity analysis; not selected by OPES hierarchy"
        ),
      TRUE ~
        "Excluded from kappa pooling: only non-kappa metric(s) reported"
    )
  ) %>%
  arrange(
    factor(analytic_stratum, levels = primary_strata),
    year,
    first_author
  ) %>%
  select(
    analytic_stratum, study_id, first_author, year, publication_type,
    retained_comparisons, retained_metric_types, eligible_kappa_effects,
    selected_opes_effects, primary_model_status
  )
write_csv(
  study_stratum_accounting,
  file.path(out_dir, "per_stratum_study_accounting.csv")
)

version_audit <- stud %>%
  filter(study_id %in% included_study_ids) %>%
  left_join(
    active_comp %>%
      group_by(study_id) %>%
      summarise(
        analytic_standardized_versions = paste(
          sort(unique(na.omit(c(std_version_1, std_version_2)))),
          collapse = "; "
        ),
        .groups = "drop"
      ),
    by = "study_id"
  ) %>%
  mutate(
    submitted_table_display_error = study_id %in% c(
      "Rosenkrantz_2014", "Seppala_2014", "Rocca_2016",
      "Ragel_2016", "Sanz_2016", "Pitra_2018",
      "Shaish_2019", "Lerchbaumer_2020", "Lucocq_2021"
    ),
    audit_disposition = if_else(
      submitted_table_display_error,
      paste(
        "Corrected Table 1 display from v2019 to original/pre-2019;",
        "analytic standardized coding was already original"
      ),
      "No version-display correction required"
    )
  ) %>%
  select(
    study_id, first_author, year, bosniak_version,
    analytic_standardized_versions, submitted_table_display_error,
    audit_disposition
  ) %>%
  arrange(year, first_author)
write_csv(version_audit, file.path(out_dir, "version_audit_79_studies.csv"))

cat("\nRevision sensitivity-analysis summary\n")
print(analysis_summary, n = Inf, width = Inf)
cat("\nVariance-source counts\n")
print(variance_counts, n = Inf, width = Inf)
cat("\nWeight-scheme counts\n")
print(weight_counts, n = Inf, width = Inf)


## ------------------------------------------------------------------
## Version analyses
## ------------------------------------------------------------------

version_ct <- ir %>%
  filter(
    std_modality_1 == "CT",
    std_version_1 %in% c("original", "v2019")
  ) %>%
  mutate(is_v2019 = if_else(std_version_1 == "v2019", 1, 0))

version_ct_model <- rma.mv(
  yi_z ~ is_v2019,
  vi_z,
  random = list(~ 1 | dep_id, ~ 1 | comparison_id),
  data = version_ct,
  method = "REML"
)
version_ct_robust <- coef_test(
  version_ct_model,
  vcov = "CR2",
  cluster = version_ct$dep_id
)
version_effect <- version_ct_robust[2, ]
version_crit <- safe_tcrit(version_effect$df_Satt)
version_ci_z <- version_effect$beta +
  c(-1, 1) * version_crit * version_effect$SE
version_baseline_z <- version_ct_robust$beta[1]
version_kappa_difference <- tanh(version_baseline_z + version_effect$beta) -
  tanh(version_baseline_z)
version_kappa_difference_ci <- tanh(version_baseline_z + version_ci_z) -
  tanh(version_baseline_z)

version_ct_summary <- tibble(
  analysis = "CT-only inter-reader version meta-regression",
  effects = nrow(version_ct),
  studies = n_distinct(version_ct$study_id),
  clusters = n_distinct(version_ct$dep_id),
  original_effects = sum(version_ct$is_v2019 == 0),
  v2019_effects = sum(version_ct$is_v2019 == 1),
  beta_z = version_effect$beta,
  se_z = version_effect$SE,
  df = version_effect$df_Satt,
  p_value = version_effect$p_Satt,
  ci_lower_z = version_ci_z[1],
  ci_upper_z = version_ci_z[2],
  original_reference_kappa = tanh(version_baseline_z),
  approximate_kappa_difference = version_kappa_difference,
  approximate_kappa_difference_ci_lower = version_kappa_difference_ci[1],
  approximate_kappa_difference_ci_upper = version_kappa_difference_ci[2]
)
write_csv(version_ct_summary, file.path(out_dir, "version_ct_model.csv"))

struct_order <- c(
  pooled = 1,
  single_pair = 2,
  single = 3,
  pairwise = 4,
  subgroup = 5,
  vs_reference = 6
)

version_pair_candidates <- ir %>%
  filter(std_version_1 %in% c("original", "v2019")) %>%
  group_by(study_id) %>%
  summarise(
    has_original = any(std_version_1 == "original"),
    has_v2019 = any(std_version_1 == "v2019"),
    .groups = "drop"
  ) %>%
  filter(has_original, has_v2019) %>%
  pull(study_id)

matched_pairs <- list()
unmatched_version_studies <- character()
for (sid in sort(version_pair_candidates)) {
  study_data <- ir %>%
    filter(
      study_id == sid,
      std_version_1 %in% c("original", "v2019")
    )
  original_data <- study_data %>% filter(std_version_1 == "original")
  v2019_data <- study_data %>% filter(std_version_1 == "v2019")

  common_cells <- inner_join(
    original_data %>% distinct(std_modality_1, reader_structure),
    v2019_data %>% distinct(std_modality_1, reader_structure),
    by = c("std_modality_1", "reader_structure")
  )
  if (nrow(common_cells) == 0) {
    unmatched_version_studies <- c(unmatched_version_studies, sid)
    next
  }

  common_cells <- common_cells %>%
    mutate(
      structure_rank = if_else(
        reader_structure %in% names(struct_order),
        unname(struct_order[reader_structure]),
        99
      )
    )
  cell_candidates <- vector("list", nrow(common_cells))
  for (idx in seq_len(nrow(common_cells))) {
    cell_modality <- common_cells$std_modality_1[idx]
    cell_structure <- common_cells$reader_structure[idx]
    cell_rows <- study_data %>%
      filter(
        std_modality_1 == cell_modality,
        reader_structure == cell_structure
      )
    cell_candidates[[idx]] <- tibble(
      std_modality_1 = cell_modality,
      reader_structure = cell_structure,
      structure_rank = common_cells$structure_rank[idx],
      maximum_n = max(cell_rows$nc, na.rm = TRUE)
    )
  }
  best_cell <- bind_rows(cell_candidates) %>%
    arrange(structure_rank, desc(maximum_n)) %>%
    slice(1)

  pick_one <- function(data, modality_value, structure_value) {
    data %>%
      filter(
        std_modality_1 == modality_value,
        reader_structure == structure_value
      ) %>%
      mutate(n_for_order = coalesce(nc, 0)) %>%
      arrange(desc(n_for_order), comparison_id) %>%
      slice(1)
  }

  original_row <- pick_one(
    original_data,
    best_cell$std_modality_1,
    best_cell$reader_structure
  )
  v2019_row <- pick_one(
    v2019_data,
    best_cell$std_modality_1,
    best_cell$reader_structure
  )

  matched_pairs[[length(matched_pairs) + 1]] <- tibble(
    study_id = sid,
    modality = best_cell$std_modality_1,
    reader_structure = best_cell$reader_structure,
    original_comparison_id = original_row$comparison_id,
    v2019_comparison_id = v2019_row$comparison_id,
    original_kappa = original_row$kappa_num,
    v2019_kappa = v2019_row$kappa_num,
    original_yi = original_row$yi_z,
    v2019_yi = v2019_row$yi_z,
    original_vi = original_row$vi_z,
    v2019_vi = v2019_row$vi_z,
    raw_kappa_difference = v2019_row$kappa_num - original_row$kappa_num
  )
}

matched_pairs <- bind_rows(matched_pairs) %>%
  mutate(
    yi_difference = v2019_yi - original_yi,
    vi_difference = original_vi + v2019_vi
  )
write_csv(matched_pairs, file.path(out_dir, "version_matched_pairs.csv"))

original_pooled <- rma(original_yi, original_vi, data = matched_pairs, method = "REML")
v2019_pooled <- rma(v2019_yi, v2019_vi, data = matched_pairs, method = "REML")
matched_difference_model <- rma(
  yi_difference,
  vi_difference,
  data = matched_pairs,
  method = "REML"
)
matched_original_kappa <- tanh(as.numeric(original_pooled$b[1]))
matched_v2019_kappa <- tanh(as.numeric(v2019_pooled$b[1]))
matched_difference_kappa <- tanh(
  as.numeric(original_pooled$b[1]) + as.numeric(matched_difference_model$b[1])
) - matched_original_kappa
matched_difference_kappa_ci <- tanh(
  as.numeric(original_pooled$b[1]) +
    c(matched_difference_model$ci.lb, matched_difference_model$ci.ub)
) - matched_original_kappa

matched_version_summary <- tibble(
  analysis = "Exact within-study matched version comparison",
  candidate_studies = length(version_pair_candidates),
  matched_studies = nrow(matched_pairs),
  unmatched_studies = paste(unmatched_version_studies, collapse = "; "),
  pooled_original_kappa = matched_original_kappa,
  pooled_original_ci_lower = tanh(original_pooled$ci.lb),
  pooled_original_ci_upper = tanh(original_pooled$ci.ub),
  pooled_v2019_kappa = matched_v2019_kappa,
  pooled_v2019_ci_lower = tanh(v2019_pooled$ci.lb),
  pooled_v2019_ci_upper = tanh(v2019_pooled$ci.ub),
  beta_difference_z = as.numeric(matched_difference_model$b[1]),
  beta_difference_ci_lower_z = matched_difference_model$ci.lb,
  beta_difference_ci_upper_z = matched_difference_model$ci.ub,
  p_value = matched_difference_model$pval,
  approximate_kappa_difference = matched_difference_kappa,
  approximate_kappa_difference_ci_lower = matched_difference_kappa_ci[1],
  approximate_kappa_difference_ci_upper = matched_difference_kappa_ci[2],
  raw_mean_kappa_difference = mean(matched_pairs$raw_kappa_difference)
)
write_csv(
  matched_version_summary,
  file.path(out_dir, "version_matched_summary.csv")
)


## ------------------------------------------------------------------
## Meta-regression table with confidence intervals, FDR, and counts
## ------------------------------------------------------------------

level_count_text <- function(data, variable, label_map = NULL) {
  values <- as.character(data[[variable]])
  count_data <- data %>%
    mutate(.level = values) %>%
    filter(!is.na(.level)) %>%
    group_by(.level) %>%
    summarise(
      effects = n(),
      studies = n_distinct(study_id),
      clusters = n_distinct(dep_id),
      .groups = "drop"
    )
  if (!is.null(label_map)) {
    count_data$.level <- if_else(
      count_data$.level %in% names(label_map),
      unname(label_map[count_data$.level]),
      count_data$.level
    )
  }
  paste(
    sprintf(
      "%s: %d effects/%d studies/%d clusters",
      count_data$.level,
      count_data$effects,
      count_data$studies,
      count_data$clusters
    ),
    collapse = "; "
  )
}

run_meta_regression <- function(
    data,
    formula,
    stratum,
    moderator,
    level_variable = NULL,
    level_label_map = NULL,
    contrast_label_map = NULL) {
  if (nrow(data) < 3 || n_distinct(data$dep_id) < 3) return(NULL)
  model <- rma.mv(
    formula,
    V = data$vi_z,
    random = list(~ 1 | dep_id, ~ 1 | comparison_id),
    data = data,
    method = "REML"
  )
  robust <- coef_test(model, vcov = "CR2", cluster = data$dep_id)
  result <- tibble(
    term = rownames(robust),
    beta = robust$beta,
    se = robust$SE,
    df = robust$df_Satt,
    p_raw = robust$p_Satt
  ) %>%
    filter(term != "intrcpt") %>%
    mutate(
      critical_value = safe_tcrit(df),
      ci_lower = beta - critical_value * se,
      ci_upper = beta + critical_value * se,
      stratum = stratum,
      moderator = moderator,
      contrast = term,
      effects = nrow(data),
      studies = n_distinct(data$study_id),
      clusters = n_distinct(data$dep_id),
      level_counts = if (is.null(level_variable)) {
        sprintf(
          "Overall: %d effects/%d studies/%d clusters",
          nrow(data),
          n_distinct(data$study_id),
          n_distinct(data$dep_id)
        )
      } else {
        level_count_text(data, level_variable, level_label_map)
      }
    )
  if (!is.null(contrast_label_map)) {
    result$contrast <- if_else(
      result$term %in% names(contrast_label_map),
      unname(contrast_label_map[result$term]),
      result$term
    )
  }
  result %>%
    select(
      stratum, moderator, contrast, effects, studies, clusters,
      level_counts, beta, se, ci_lower, ci_upper, df, p_raw
    )
}

meta_parts <- list()

ir_ncat <- ir %>%
  filter(ncat_group %in% c("2-3", "4-5")) %>%
  mutate(is_full = if_else(ncat_group == "4-5", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_ncat, yi_z ~ is_full, "Inter-reader", "No. of categories",
  "is_full", c("0" = "2-3", "1" = "4-5"),
  c(is_full = "4-5 vs 2-3")
)))

ir_weight <- ir %>%
  filter(wt_group %in% c("weighted", "unweighted")) %>%
  mutate(is_weighted = if_else(wt_group == "weighted", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_weight, yi_z ~ is_weighted, "Inter-reader", "Weight scheme",
  "is_weighted", c("0" = "Unweighted", "1" = "Weighted"),
  c(is_weighted = "Weighted vs unweighted")
)))

ir_modality <- ir %>%
  filter(mod_group %in% c("CT", "MRI", "CEUS_US")) %>%
  mutate(modality = factor(mod_group, levels = c("CT", "MRI", "CEUS_US")))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_modality, yi_z ~ modality, "Inter-reader", "Imaging modality",
  "mod_group", NULL,
  c(
    modalityMRI = "MRI vs CT",
    modalityCEUS_US = "CEUS/US vs CT"
  )
)))

meta_parts <- append(meta_parts, list(run_meta_regression(
  version_ct, yi_z ~ is_v2019, "Inter-reader", "Bosniak version (CT-only)",
  "is_v2019", c("0" = "Original", "1" = "v2019"),
  c(is_v2019 = "v2019 vs original")
)))

ir_publication <- ir %>%
  filter(pub_group %in% c("fulltext", "abstract")) %>%
  mutate(is_abstract = if_else(pub_group == "abstract", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_publication, yi_z ~ is_abstract, "Inter-reader", "Publication type",
  "is_abstract", c("0" = "Full text", "1" = "Abstract"),
  c(is_abstract = "Abstract vs full text")
)))

meta_parts <- append(meta_parts, list(run_meta_regression(
  ir %>% filter(!is.na(log_nc)),
  yi_z ~ log_nc, "Inter-reader", "log(sample size)"
)))

ir_readers <- ir %>%
  mutate(
    reader_count_group = case_when(
      n_readers_num == 2 ~ "2",
      n_readers_num >= 3 ~ "3plus",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(reader_count_group)) %>%
  mutate(is_3plus = if_else(reader_count_group == "3plus", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_readers, yi_z ~ is_3plus, "Inter-reader", "No. of readers",
  "is_3plus", c("0" = "2", "1" = "3 or more"),
  c(is_3plus = "3 or more vs 2")
)))

meta_parts <- append(meta_parts, list(run_meta_regression(
  ir %>% filter(!is.na(exp_mean_num)),
  yi_z ~ exp_mean_num, "Inter-reader", "Mean reader experience"
)))

ir_specialty <- ir %>%
  filter(!is.na(spec_group)) %>%
  mutate(is_subspecialty = if_else(spec_group == "subspecialty", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_specialty, yi_z ~ is_subspecialty,
  "Inter-reader", "Reader subspecialization",
  "is_subspecialty",
  c("0" = "Other", "1" = "Subspecialty"),
  c(is_subspecialty = "Subspecialty vs other")
)))

ir_blinding <- ir %>%
  filter(!is.na(blind_group)) %>%
  mutate(blinding_reported_binary = if_else(blind_group == "reported", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_blinding, yi_z ~ blinding_reported_binary,
  "Inter-reader", "Blinding reported",
  "blinding_reported_binary",
  c("0" = "Not reported", "1" = "Reported"),
  c(blinding_reported_binary = "Reported vs not reported")
)))

ir_structure <- ir %>%
  filter(!is.na(structure_group)) %>%
  mutate(is_pooled = if_else(structure_group == "pooled", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_structure, yi_z ~ is_pooled, "Inter-reader", "Reader structure",
  "is_pooled", c("0" = "Other", "1" = "Pooled"),
  c(is_pooled = "Pooled vs other")
)))

meta_parts <- append(meta_parts, list(run_meta_regression(
  ir %>% filter(!is.na(year_num)),
  yi_z ~ year_num, "Inter-reader", "Publication year"
)))

ir_design <- ir %>%
  filter(design_group %in% c("retrospective", "prospective")) %>%
  mutate(is_prospective = if_else(design_group == "prospective", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_design, yi_z ~ is_prospective, "Inter-reader", "Study design",
  "is_prospective",
  c("0" = "Retrospective", "1" = "Prospective"),
  c(is_prospective = "Prospective vs retrospective")
)))

ir_region <- ir %>%
  filter(region_group %in% c("North_America", "Asia", "Europe")) %>%
  mutate(
    region_factor = factor(
      region_group,
      levels = c("North_America", "Asia", "Europe")
    )
  )
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_region, yi_z ~ region_factor, "Inter-reader", "Region",
  "region_group", NULL,
  c(
    region_factorAsia = "Asia vs North America",
    region_factorEurope = "Europe vs North America"
  )
)))

ir_rep <- ir %>%
  mutate(rep_yes = if_else(representative_spectrum == "Yes", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  ir_rep, yi_z ~ rep_yes, "Inter-reader", "Representative spectrum",
  "rep_yes",
  c("0" = "No/unclear", "1" = "Yes"),
  c(rep_yes = "Yes vs no/unclear")
)))

im_ncat <- im %>%
  filter(ncat_group %in% c("2-3", "4-5")) %>%
  mutate(is_full = if_else(ncat_group == "4-5", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_ncat, yi_z ~ is_full, "Inter-modality", "No. of categories",
  "is_full", c("0" = "2-3", "1" = "4-5"),
  c(is_full = "4-5 vs 2-3")
)))

im_pair <- im %>%
  filter(mod_pair %in% c("CT_MRI", "CT_CEUS_US")) %>%
  mutate(
    pair_factor = factor(mod_pair, levels = c("CT_MRI", "CT_CEUS_US"))
  )
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_pair, yi_z ~ pair_factor, "Inter-modality", "Modality pair",
  "mod_pair", NULL,
  c(pair_factorCT_CEUS_US = "CT-US-based vs CT-MRI")
)))

meta_parts <- append(meta_parts, list(run_meta_regression(
  im %>% filter(!is.na(log_nc)),
  yi_z ~ log_nc, "Inter-modality", "log(sample size)"
)))

im_publication <- im %>%
  filter(pub_group %in% c("fulltext", "abstract")) %>%
  mutate(is_abstract = if_else(pub_group == "abstract", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_publication, yi_z ~ is_abstract, "Inter-modality", "Publication type",
  "is_abstract", c("0" = "Full text", "1" = "Abstract"),
  c(is_abstract = "Abstract vs full text")
)))

meta_parts <- append(meta_parts, list(run_meta_regression(
  im %>% filter(!is.na(exp_mean_num)),
  yi_z ~ exp_mean_num, "Inter-modality", "Mean reader experience"
)))

im_specialty <- im %>%
  filter(!is.na(spec_group)) %>%
  mutate(is_subspecialty = if_else(spec_group == "subspecialty", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_specialty, yi_z ~ is_subspecialty,
  "Inter-modality", "Reader subspecialization",
  "is_subspecialty",
  c("0" = "Other", "1" = "Subspecialty"),
  c(is_subspecialty = "Subspecialty vs other")
)))

im_blinding <- im %>%
  filter(!is.na(blind_group)) %>%
  mutate(blinding_reported_binary = if_else(blind_group == "reported", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_blinding, yi_z ~ blinding_reported_binary,
  "Inter-modality", "Blinding reported",
  "blinding_reported_binary",
  c("0" = "Not reported", "1" = "Reported"),
  c(blinding_reported_binary = "Reported vs not reported")
)))

meta_parts <- append(meta_parts, list(run_meta_regression(
  im %>% filter(!is.na(year_num)),
  yi_z ~ year_num, "Inter-modality", "Publication year"
)))

im_weight <- im %>%
  filter(wt_group %in% c("weighted", "unweighted")) %>%
  mutate(is_weighted = if_else(wt_group == "weighted", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_weight, yi_z ~ is_weighted, "Inter-modality", "Weight scheme",
  "is_weighted", c("0" = "Unweighted", "1" = "Weighted"),
  c(is_weighted = "Weighted vs unweighted")
)))

im_readers <- im %>%
  mutate(
    reader_count_group = case_when(
      n_readers_num == 1 ~ "1",
      n_readers_num == 2 ~ "2",
      n_readers_num >= 3 ~ "3plus",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(reader_count_group)) %>%
  mutate(
    reader_count_factor = factor(
      reader_count_group,
      levels = c("1", "2", "3plus")
    )
  )
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_readers, yi_z ~ reader_count_factor,
  "Inter-modality", "No. of readers",
  "reader_count_group", NULL,
  c(
    reader_count_factor2 = "2 vs 1",
    reader_count_factor3plus = "3 or more vs 1"
  )
)))

im_design <- im %>%
  filter(design_group %in% c("retrospective", "prospective")) %>%
  mutate(is_prospective = if_else(design_group == "prospective", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_design, yi_z ~ is_prospective, "Inter-modality", "Study design",
  "is_prospective",
  c("0" = "Retrospective", "1" = "Prospective"),
  c(is_prospective = "Prospective vs retrospective")
)))

im_region <- im %>%
  filter(region_group %in% c("North_America", "Asia", "Europe")) %>%
  mutate(
    region_factor = factor(
      region_group,
      levels = c("North_America", "Asia", "Europe")
    )
  )
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_region, yi_z ~ region_factor, "Inter-modality", "Region",
  "region_group", NULL,
  c(
    region_factorAsia = "Asia vs North America",
    region_factorEurope = "Europe vs North America"
  )
)))

im_rep <- im %>%
  mutate(rep_yes = if_else(representative_spectrum == "Yes", 1, 0))
meta_parts <- append(meta_parts, list(run_meta_regression(
  im_rep, yi_z ~ rep_yes, "Inter-modality", "Representative spectrum",
  "rep_yes",
  c("0" = "No/unclear", "1" = "Yes"),
  c(rep_yes = "Yes vs no/unclear")
)))

meta_regression_results <- bind_rows(meta_parts) %>%
  mutate(
    p_fdr_bh = p.adjust(p_raw, method = "BH"),
    df_interpretation = case_when(
      stratum == "Inter-modality" &
        moderator == "Publication type" &
        grepl("Abstract: 1 effects/1 studies/1 clusters", level_counts) ~
        "Non-interpretable (singleton moderator level)",
      df < 4 ~ "Non-interpretable (df<4)",
      df < 10 ~ "Fragile (df 4-<10)",
      TRUE ~ "Exploratory"
    )
  )
write_csv(
  meta_regression_results,
  file.path(out_dir, "meta_regression_revision.csv")
)

cat("\nVersion CT-only model\n")
print(version_ct_summary, n = Inf, width = Inf)
cat("\nMatched version comparison\n")
print(matched_version_summary, n = Inf, width = Inf)
cat("\nMeta-regression results with BH-FDR\n")
print(meta_regression_results, n = Inf, width = Inf)

## Leave-one-study-out checks for every moderator contrast meeting the
## protocol's nominal p<0.05 threshold. The full models are first checked
## against meta_regression_results so that the deletion analyses cannot drift
## from the reported model definitions.

fit_loso_target <- function(data, formula, target_term) {
  captured_warnings <- character()
  tryCatch(
    withCallingHandlers({
      if (nrow(data) < 3 || n_distinct(data$dep_id) < 3) {
        stop("fewer than three effects or dependency clusters")
      }
      model <- rma.mv(
        formula,
        V = data$vi_z,
        random = list(~ 1 | dep_id, ~ 1 | comparison_id),
        data = data,
        method = "REML"
      )
      robust <- coef_test(model, vcov = "CR2", cluster = data$dep_id)
      robust_df <- as.data.frame(robust)
      robust_df$term <- rownames(robust_df)
      target <- robust_df[robust_df$term == target_term, , drop = FALSE]
      if (nrow(target) != 1) {
        stop(sprintf("target term '%s' was not estimable", target_term))
      }
      if (!all(is.finite(c(
        target$beta, target$SE, target$df_Satt, target$p_Satt
      )))) {
        stop(sprintf("target term '%s' had non-finite CR2 output", target_term))
      }
      list(
        ok = TRUE,
        beta = as.numeric(target$beta),
        se = as.numeric(target$SE),
        df = as.numeric(target$df_Satt),
        p_raw = as.numeric(target$p_Satt),
        warnings = paste(unique(captured_warnings), collapse = " | "),
        error = ""
      )
    }, warning = function(w) {
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      list(
        ok = FALSE,
        beta = NA_real_,
        se = NA_real_,
        df = NA_real_,
        p_raw = NA_real_,
        warnings = paste(unique(captured_warnings), collapse = " | "),
        error = conditionMessage(e)
      )
    }
  )
}

loso_specs <- list(
  list(
    analysis_id = "IR_blinding_reported",
    stratum = "Inter-reader",
    moderator = "Blinding reported",
    contrast = "Reported vs not reported",
    data = ir_blinding,
    formula = yi_z ~ blinding_reported_binary,
    target_term = "blinding_reported_binary",
    target_level_studies = NA_character_
  ),
  list(
    analysis_id = "IR_region_Europe_vs_North_America",
    stratum = "Inter-reader",
    moderator = "Region",
    contrast = "Europe vs North America",
    data = ir_region,
    formula = yi_z ~ region_factor,
    target_term = "region_factorEurope",
    target_level_studies = NA_character_
  ),
  list(
    analysis_id = "IM_categories_4_5_vs_2_3",
    stratum = "Inter-modality",
    moderator = "No. of categories",
    contrast = "4-5 vs 2-3",
    data = im_ncat,
    formula = yi_z ~ is_full,
    target_term = "is_full",
    target_level_studies = NA_character_
  ),
  list(
    analysis_id = "IM_modality_pair_CT_US_vs_CT_MRI",
    stratum = "Inter-modality",
    moderator = "Modality pair",
    contrast = "CT-US-based vs CT-MRI",
    data = im_pair,
    formula = yi_z ~ pair_factor,
    target_term = "pair_factorCT_CEUS_US",
    target_level_studies = NA_character_
  ),
  list(
    analysis_id = "IM_log_sample_size",
    stratum = "Inter-modality",
    moderator = "log(sample size)",
    contrast = "log_nc",
    data = im %>% filter(!is.na(log_nc)),
    formula = yi_z ~ log_nc,
    target_term = "log_nc",
    target_level_studies = NA_character_
  ),
  list(
    analysis_id = "IM_publication_type_abstract_vs_fulltext",
    stratum = "Inter-modality",
    moderator = "Publication type",
    contrast = "Abstract vs full text",
    data = im_publication,
    formula = yi_z ~ is_abstract,
    target_term = "is_abstract",
    target_level_studies = paste(
      sort(unique(im_publication$study_id[im_publication$is_abstract == 1])),
      collapse = ";"
    )
  ),
  list(
    analysis_id = "IM_readers_3plus_vs_1",
    stratum = "Inter-modality",
    moderator = "No. of readers",
    contrast = "3 or more vs 1",
    data = im_readers,
    formula = yi_z ~ reader_count_factor,
    target_term = "reader_count_factor3plus",
    target_level_studies = NA_character_
  ),
  list(
    analysis_id = "IM_representative_spectrum",
    stratum = "Inter-modality",
    moderator = "Representative spectrum",
    contrast = "Yes vs no/unclear",
    data = im_rep,
    formula = yi_z ~ rep_yes,
    target_term = "rep_yes",
    target_level_studies = NA_character_
  )
)

canonical_nominal <- meta_regression_results %>%
  filter(p_raw < 0.05) %>%
  select(stratum, moderator, contrast, beta, se, df, p_raw, df_interpretation)

if (nrow(canonical_nominal) != 8) {
  stop(sprintf(
    "Expected eight moderator contrasts with nominal p<0.05, found %d.",
    nrow(canonical_nominal)
  ))
}

spec_keys <- bind_rows(lapply(loso_specs, function(x) {
  tibble(
    stratum = x$stratum,
    moderator = x$moderator,
    contrast = x$contrast
  )
}))
key_columns <- c("stratum", "moderator", "contrast")
if (
  nrow(anti_join(canonical_nominal, spec_keys, by = key_columns)) > 0 ||
    nrow(anti_join(spec_keys, canonical_nominal, by = key_columns)) > 0
) {
  stop("LOSO specifications do not match the nominal-p meta-regression rows.")
}

loso_validation_parts <- list()
loso_iteration_parts <- list()

for (spec in loso_specs) {
  full_fit <- fit_loso_target(spec$data, spec$formula, spec$target_term)
  if (!isTRUE(full_fit$ok)) {
    stop(sprintf(
      "Full model failed for %s: %s",
      spec$analysis_id,
      full_fit$error
    ))
  }

  canonical_row <- canonical_nominal %>%
    filter(
      stratum == spec$stratum,
      moderator == spec$moderator,
      contrast == spec$contrast
    )
  if (nrow(canonical_row) != 1) {
    stop(sprintf("Canonical row lookup failed for %s.", spec$analysis_id))
  }

  loso_validation_parts[[length(loso_validation_parts) + 1]] <- tibble(
    analysis_id = spec$analysis_id,
    stratum = spec$stratum,
    moderator = spec$moderator,
    contrast = spec$contrast,
    target_term = spec$target_term,
    effects = nrow(spec$data),
    studies = n_distinct(spec$data$study_id),
    clusters = n_distinct(spec$data$dep_id),
    beta_recomputed = full_fit$beta,
    beta_canonical = canonical_row$beta,
    beta_abs_diff = abs(full_fit$beta - canonical_row$beta),
    se_recomputed = full_fit$se,
    se_canonical = canonical_row$se,
    se_abs_diff = abs(full_fit$se - canonical_row$se),
    df_recomputed = full_fit$df,
    df_canonical = canonical_row$df,
    df_abs_diff = abs(full_fit$df - canonical_row$df),
    p_recomputed = full_fit$p_raw,
    p_canonical = canonical_row$p_raw,
    p_abs_diff = abs(full_fit$p_raw - canonical_row$p_raw),
    canonical_interpretation = canonical_row$df_interpretation,
    target_level_studies = spec$target_level_studies,
    full_model_warnings = full_fit$warnings
  )

  for (omitted_study in sort(unique(spec$data$study_id))) {
    deletion_data <- spec$data %>% filter(study_id != omitted_study)
    deletion_fit <- fit_loso_target(
      deletion_data,
      spec$formula,
      spec$target_term
    )
    loso_iteration_parts[[length(loso_iteration_parts) + 1]] <- tibble(
      analysis_id = spec$analysis_id,
      stratum = spec$stratum,
      moderator = spec$moderator,
      contrast = spec$contrast,
      omitted_study = omitted_study,
      removed_effects = sum(spec$data$study_id == omitted_study),
      remaining_effects = nrow(deletion_data),
      remaining_studies = n_distinct(deletion_data$study_id),
      remaining_clusters = n_distinct(deletion_data$dep_id),
      successful = deletion_fit$ok,
      beta = deletion_fit$beta,
      se = deletion_fit$se,
      df = deletion_fit$df,
      p_raw = deletion_fit$p_raw,
      nominal_p_lt_0_05 = ifelse(
        deletion_fit$ok,
        deletion_fit$p_raw < 0.05,
        NA
      ),
      direction_consistent = ifelse(
        deletion_fit$ok,
        sign(deletion_fit$beta) == sign(full_fit$beta),
        NA
      ),
      warnings = deletion_fit$warnings,
      failure_reason = deletion_fit$error
    )
  }
}

loso_validation <- bind_rows(loso_validation_parts)
loso_iterations <- bind_rows(loso_iteration_parts)

validation_tolerance <- 1e-10
if (any(
  loso_validation$beta_abs_diff > validation_tolerance |
    loso_validation$se_abs_diff > validation_tolerance |
    loso_validation$df_abs_diff > validation_tolerance |
    loso_validation$p_abs_diff > validation_tolerance
)) {
  stop("At least one LOSO full model did not match the canonical result.")
}

loso_summary <- loso_iterations %>%
  left_join(
    loso_validation %>%
      select(
        analysis_id,
        full_beta = beta_recomputed,
        full_p_raw = p_recomputed,
        canonical_interpretation,
        target_level_studies
      ),
    by = "analysis_id"
  ) %>%
  group_by(
    analysis_id, stratum, moderator, contrast,
    full_beta, full_p_raw, canonical_interpretation,
    target_level_studies
  ) %>%
  summarise(
    planned_iterations = n(),
    successful_iterations = sum(successful),
    non_estimable_iterations = sum(!successful),
    direction_consistent_n = sum(direction_consistent, na.rm = TRUE),
    direction_consistent_denominator = sum(successful),
    nominal_p_lt_0_05_n = sum(nominal_p_lt_0_05, na.rm = TRUE),
    nominal_p_lt_0_05_denominator = sum(successful),
    beta_min = ifelse(any(successful), min(beta[successful]), NA_real_),
    beta_max = ifelse(any(successful), max(beta[successful]), NA_real_),
    p_raw_min = ifelse(any(successful), min(p_raw[successful]), NA_real_),
    p_raw_max = ifelse(any(successful), max(p_raw[successful]), NA_real_),
    non_estimable_studies = paste(omitted_study[!successful], collapse = ";"),
    failure_reasons = paste(unique(failure_reason[!successful]), collapse = " | "),
    warning_iterations = sum(warnings != ""),
    .groups = "drop"
  )

write_csv(
  loso_iterations,
  file.path(out_dir, "loso_nominal_moderator_iterations.csv")
)
write_csv(
  loso_summary,
  file.path(out_dir, "loso_nominal_moderator_summary.csv")
)
write_csv(
  loso_validation,
  file.path(out_dir, "loso_nominal_moderator_validation.csv")
)

## Focused leave-one-abstract-study-out checks for the inter-reader
## publication-type contrast. These are separate from the nominal-p set above
## because the full inter-reader publication-type contrast had p=0.272.

publication_check_parts <- list()
publication_check_specs <- c("Full model", sort(unique(
  ir_publication$study_id[ir_publication$is_abstract == 1]
)))

for (check_label in publication_check_specs) {
  omitted_study <- if (check_label == "Full model") NA_character_ else check_label
  check_data <- if (is.na(omitted_study)) {
    ir_publication
  } else {
    ir_publication %>% filter(study_id != omitted_study)
  }
  check_fit <- fit_loso_target(
    check_data,
    yi_z ~ is_abstract,
    "is_abstract"
  )
  abstract_data <- check_data %>% filter(is_abstract == 1)
  abstract_studies <- n_distinct(abstract_data$study_id)
  abstract_clusters <- n_distinct(abstract_data$dep_id)
  publication_check_parts[[length(publication_check_parts) + 1]] <- tibble(
    analysis = check_label,
    omitted_study = omitted_study,
    effects = nrow(check_data),
    studies = n_distinct(check_data$study_id),
    clusters = n_distinct(check_data$dep_id),
    abstract_effects = nrow(abstract_data),
    abstract_studies = abstract_studies,
    abstract_clusters = abstract_clusters,
    successful = check_fit$ok,
    beta = check_fit$beta,
    se = check_fit$se,
    df = check_fit$df,
    p_raw = check_fit$p_raw,
    interpretation = case_when(
      !check_fit$ok ~ "Non-estimable",
      abstract_studies < 2 || abstract_clusters < 2 ~
        "Non-interpretable (singleton abstract level)",
      check_fit$df < 4 ~ "Non-interpretable (df<4)",
      check_fit$df < 10 ~ "Fragile (df 4-<10)",
      TRUE ~ "Exploratory"
    ),
    warnings = check_fit$warnings,
    failure_reason = check_fit$error
  )
}

inter_reader_publication_checks <- bind_rows(publication_check_parts)
write_csv(
  inter_reader_publication_checks,
  file.path(out_dir, "inter_reader_publication_type_leave_one_abstract.csv")
)

cat("\nNominal-p moderator leave-one-study-out summary\n")
print(loso_summary, n = Inf, width = Inf)
cat("\nInter-reader publication-type leave-one-abstract checks\n")
print(inter_reader_publication_checks, n = Inf, width = Inf)
