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
 ÷o6¶‰žËkºwµç@¹}É•…‘•ÉÍ}¹Õ´€ôô€Äø€ˆÄˆ°(€€€€€¹}É•…‘•ÉÍ}¹Õ´€ôô€Èø€ˆÈˆ°(€€€€€¹}É•…‘•ÉÍ}¹Õ´€øô€Ìø€ˆÍÁ±ÕÌˆ°(€€€€€QIUø9}¡…É…Ñ•É|(€€€€¤(€€¤€”ø”(€™¥±Ñ•È …¥Ì¹¹„¡É•…‘•É}½Õ¹Ñ}É½ÕÀ¤¤€”ø”(€µÕÑ…Ñ” (€€€É•…‘•É}½Õ¹Ñ}™…Ñ½È€ô™…Ñ½È (€€€€€É•…‘•É}½Õ¹Ñ}É½ÕÀ°(€€€€€±•Ù•±Ì€ôŒ ˆÄˆ°€ˆÈˆ°€ˆÍÁ±ÕÌˆ¤(€€€€¤(€€¤)µ•Ñ…}Á…ÉÑÌ€ð´…ÁÁ•¹¡µ•Ñ…}Á…ÉÑÌ°±¥ÍÐ¡ÉÕ¹}µ•Ñ…}É•É•ÍÍ¥½¸ (€¥µ}É•…‘•ÉÌ°å¥}èøÉ•…‘•É}½Õ¹Ñ}™…Ñ½È°(€€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°€‰9¼¸½˜É•…‘•ÉÌˆ°(€€‰É•…‘•É}½Õ¹Ñ}É½ÕÀˆ°9U10°(€Œ (€€€É•…‘•É}½Õ¹Ñ}™…Ñ½ÈÈ€ô€ˆÈÙÌ€Äˆ°(€€€É•…‘•É}½Õ¹Ñ}™…Ñ½ÈÍÁ±ÕÌ€ô€ˆÌ½Èµ½É”ÙÌ€Äˆ(€€¤(¤¤¤()¥µ}‘•Í¥¸€ð´¥´€”ø”(€™¥±Ñ•È¡‘•Í¥¹}É½ÕÀ€•¥¸”Œ ‰É•ÑÉ½ÍÁ•Ñ¥Ù”ˆ°€‰ÁÉ½ÍÁ•Ñ¥Ù”ˆ¤¤€”ø”(€µÕÑ…Ñ”¡¥Í}ÁÉ½ÍÁ•Ñ¥Ù”€ô¥™}•±Í”¡‘•Í¥¹}É½ÕÀ€ôô€‰ÁÉ½ÍÁ•Ñ¥Ù”ˆ°€Ä°€À¤¤)µ•Ñ…}Á…ÉÑÌ€ð´…ÁÁ•¹¡µ•Ñ…}Á…ÉÑÌ°±¥ÍÐ¡ÉÕ¹}µ•Ñ…}É•É•ÍÍ¥½¸ (€¥µ}‘•Í¥¸°å¥}èø¥Í}ÁÉ½ÍÁ•Ñ¥Ù”°€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°€‰MÑÕ‘ä‘•Í¥¸ˆ°(€€‰¥Í}ÁÉ½ÍÁ•Ñ¥Ù”ˆ°(€Œ ˆÀˆ€ô€‰I•ÑÉ½ÍÁ•Ñ¥Ù”ˆ°€ˆÄˆ€ô€‰AÉ½ÍÁ•Ñ¥Ù”ˆ¤°(€Œ¡¥Í}ÁÉ½ÍÁ•Ñ¥Ù”€ô€‰AÉ½ÍÁ•Ñ¥Ù”ÙÌÉ•ÑÉ½ÍÁ•Ñ¥Ù”ˆ¤(¤¤¤()¥µ}É•¥½¸€ð´¥´€”ø”(€™¥±Ñ•È¡É•¥½¹}É½ÕÀ€•¥¸”Œ ‰9½ÉÑ¡}µ•É¥„ˆ°€‰Í¥„ˆ°€‰ÕÉ½Á”ˆ¤¤€”ø”(€µÕÑ…Ñ” (€€€É•¥½¹}™…Ñ½È€ô™…Ñ½È (€€€€€É•¥½¹}É½ÕÀ°(€€€€€±•Ù•±Ì€ôŒ ‰9½ÉÑ¡}µ•É¥„ˆ°€‰Í¥„ˆ°€‰ÕÉ½Á”ˆ¤(€€€€¤(€€¤)µ•Ñ…}Á…ÉÑÌ€ð´…ÁÁ•¹¡µ•Ñ…}Á…ÉÑÌ°±¥ÍÐ¡ÉÕ¹}µ•Ñ…}É•É•ÍÍ¥½¸ (€¥µ}É•¥½¸°å¥}èøÉ•¥½¹}™…Ñ½È°€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°€‰I•¥½¸ˆ°(€€‰É•¥½¹}É½ÕÀˆ°9U10°(€Œ (€€€É•¥½¹}™…Ñ½ÉÍ¥„€ô€‰Í¥„ÙÌ9½ÉÑ µ•É¥„ˆ°(€€€É•¥½¹}™…Ñ½ÉÕÉ½Á”€ô€‰ÕÉ½Á”ÙÌ9½ÉÑ µ•É¥„ˆ(€€¤(¤¤¤()¥µ}É•À€ð´¥´€”ø”(€µÕÑ…Ñ”¡É•Á}å•Ì€ô¥™}•±Í”¡É•ÁÉ•Í•¹Ñ…Ñ¥Ù•}ÍÁ•ÑÉÕ´€ôô€‰e•Ìˆ°€Ä°€À¤¤)µ•Ñ…}Á…ÉÑÌ€ð´…ÁÁ•¹¡µ•Ñ…}Á…ÉÑÌ°±¥ÍÐ¡ÉÕ¹}µ•Ñ…}É•É•ÍÍ¥½¸ (€¥µ}É•À°å¥}èøÉ•Á}å•Ì°€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°€‰I•ÁÉ•Í•¹Ñ…Ñ¥Ù”ÍÁ•ÑÉÕ´ˆ°(€€‰É•Á}å•Ìˆ°(€Œ ˆÀˆ€ô€‰9¼½Õ¹±•…Èˆ°€ˆÄˆ€ô€‰e•Ìˆ¤°(€Œ¡É•Á}å•Ì€ô€‰e•ÌÙÌ¹¼½Õ¹±•…Èˆ¤(¤¤¤()µ•Ñ…}É•É•ÍÍ¥½¹}É•ÍÕ±ÑÌ€ð´‰¥¹‘}É½ÝÌ¡µ•Ñ…}Á…ÉÑÌ¤€”ø”(€µÕÑ…Ñ” (€€€Á}™‘É}‰ €ôÀ¹…‘©ÕÍÐ¡Á}É…Ü°µ•Ñ¡½€ô€‰	 ˆ¤°(€€€‘™}¥¹Ñ•ÉÁÉ•Ñ…Ñ¥½¸€ô…Í•}Ý¡•¸ (€€€€€ÍÑÉ…ÑÕ´€ôô€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ€˜(€€€€€€€µ½‘•É…Ñ½È€ôô€‰AÕ‰±¥…Ñ¥½¸ÑåÁ”ˆ€˜(€€€€€€€É•Á° ‰‰ÍÑÉ…Ðè€Ä•™™•ÑÌ¼ÄÍÑÕ‘¥•Ì¼Ä±ÕÍÑ•ÉÌˆ°±•Ù•±}½Õ¹ÑÌ¤ø(€€€€€€€€‰9½¸µ¥¹Ñ•ÉÁÉ•Ñ…‰±”€¡Í¥¹±•Ñ½¸µ½‘•É…Ñ½È±•Ù•°¤ˆ°(€€€€€‘˜€ð€Ðø€‰9½¸µ¥¹Ñ•ÉÁÉ•Ñ…‰±”€¡‘˜ðÐ¤ˆ°(€€€€€‘˜€ð€ÄÀø€‰É…¥±”€¡‘˜€Ð´ðÄÀ¤ˆ°(€€€€€QIUø€‰áÁ±½É…Ñ½Éäˆ(€€€€¤(€€¤)ÝÉ¥Ñ•}ÍØ (€µ•Ñ…}É•É•ÍÍ¥½¹}É•ÍÕ±ÑÌ°(€™¥±”¹Á…Ñ ¡½ÕÑ}‘¥È°€‰µ•Ñ…}É•É•ÍÍ¥½¹}É•Ù¥Í¥½¸¹ÍØˆ¤(¤()…Ð ‰q¹Y•ÉÍ¥½¸Pµ½¹±äµ½‘•±q¸ˆ¤)ÁÉ¥¹Ð¡Ù•ÉÍ¥½¹}Ñ}ÍÕµµ…Éä°¸€ô%¹˜°Ý¥‘Ñ €ô%¹˜¤)…Ð ‰q¹5…Ñ¡•Ù•ÉÍ¥½¸½µÁ…É¥Í½¹q¸ˆ¤)ÁÉ¥¹Ð¡µ…Ñ¡•‘}Ù•ÉÍ¥½¹}ÍÕµµ…Éä°¸€ô%¹˜°Ý¥‘Ñ €ô%¹˜¤)…Ð ‰q¹5•Ñ„µÉ•É•ÍÍ¥½¸É•ÍÕ±ÑÌÝ¥Ñ 	 µIq¸ˆ¤)ÁÉ¥¹Ð¡µ•Ñ…}É•É•ÍÍ¥½¹}É•ÍÕ±ÑÌ°¸€ô%¹˜°Ý¥‘Ñ €ô%¹˜¤((ŒŒ1•…Ù”µ½¹”µÍÑÕ‘äµ½ÕÐ¡•­Ì™½È•Ù•Éäµ½‘•É…Ñ½È½¹ÑÉ…ÍÐµ••Ñ¥¹œÑ¡”(ŒŒÁÉ½Ñ½½°Ì¹½µ¥¹…°ÀðÀ¸ÀÔÑ¡É•Í¡½±¸Q¡”™Õ±°µ½‘•±Ì…É”™¥ÉÍÐ¡•­•(ŒŒ……¥¹ÍÐµ•Ñ…}É•É•ÍÍ¥½¹}É•ÍÕ±ÑÌÍ¼Ñ¡…ÐÑ¡”‘•±•Ñ¥½¸…¹…±åÍ•Ì…¹¹½Ð‘É¥™Ð(ŒŒ™É½´Ñ¡”É•Á½ÉÑ•µ½‘•°‘•™¥¹¥Ñ¥½¹Ì¸()™¥Ñ}±½Í½}Ñ…É•Ð€ð´™Õ¹Ñ¥½¸¡‘…Ñ„°™½ÉµÕ±„°Ñ…É•Ñ}Ñ•É´¤ì(€…ÁÑÕÉ•‘}Ý…É¹¥¹Ì€ð´¡…É…Ñ•È ¤(€ÑÉå…Ñ  (€€€Ý¥Ñ¡…±±¥¹!…¹‘±•ÉÌ¡ì(€€€€€¥˜€¡¹É½Ü¡‘…Ñ„¤€ð€Ìñð¹}‘¥ÍÑ¥¹Ð¡‘…Ñ„‘‘•Á}¥¤€ð€Ì¤ì(€€€€€€€ÍÑ½À ‰™•Ý•ÈÑ¡…¸Ñ¡É•”•™™•ÑÌ½È‘•Á•¹‘•¹ä±ÕÍÑ•ÉÌˆ¤(€€€€€ô(€€€€€µ½‘•°€ð´Éµ„¹µØ (€€€€€€€™½ÉµÕ±„°(€€€€€€€X€ô‘…Ñ„‘Ù¥}è°(€€€€€€€É…¹‘½´€ô±¥ÍÐ¡ø€Äð‘•Á}¥°ø€Äð½µÁ…É¥Í½¹}¥¤°(€€€€€€€‘…Ñ„€ô‘…Ñ„°(€€€€€€€µ•Ñ¡½€ô€‰I50ˆ(€€€€€€¤(€€€€€É½‰ÕÍÐ€ð´½•™}Ñ•ÍÐ¡µ½‘•°°Ù½Ø€ô€‰HÈˆ°±ÕÍÑ•È€ô‘…Ñ„‘‘•Á}¥¤(€€€€€É½‰ÕÍÑ}‘˜€ð´…Ì¹‘…Ñ„¹™É…µ”¡É½‰ÕÍÐ¤(€€€€€É½‰ÕÍÑ}‘˜‘Ñ•É´€ð´É½Ý¹…µ•Ì¡É½‰ÕÍÑ}‘˜¤(€€€€€Ñ…É•Ð€ð´É½‰ÕÍÑ}‘™mÉ½‰ÕÍÑ}‘˜‘Ñ•É´€ôôÑ…É•Ñ}Ñ•É´°€°‘É½À€ô1Mt(€€€€€¥˜€¡¹É½Ü¡Ñ…É•Ð¤€„ô€Ä¤ì(€€€€€€€ÍÑ½À¡ÍÁÉ¥¹Ñ˜ ‰Ñ…É•ÐÑ•É´€œ•ÌœÝ…Ì¹½Ð•ÍÑ¥µ…‰±”ˆ°Ñ…É•Ñ}Ñ•É´¤¤(€€€€€ô(€€€€€¥˜€ ……±°¡¥Ì¹™¥¹¥Ñ”¡Œ (€€€€€€€Ñ…É•Ð‘‰•Ñ„°Ñ…É•Ð‘M°Ñ…É•Ð‘‘™}M…ÑÐ°Ñ…É•Ð‘Á}M…ÑÐ(€€€€€€¤¤¤¤ì(€€€€€€€ÍÑ½À¡ÍÁÉ¥¹Ñ˜ ‰Ñ…É•ÐÑ•É´€œ•Ìœ¡…¹½¸µ™¥¹¥Ñ”HÈ½ÕÑÁÕÐˆ°Ñ…É•Ñ}Ñ•É´¤¤(€€€€€ô(€€€€€±¥ÍÐ (€€€€€€€½¬€ôQIU°(€€€€€€€‰•Ñ„€ô…Ì¹¹Õµ•É¥Œ¡Ñ…É•Ð‘‰•Ñ„¤°(€€€€€€€Í”€ô…Ì¹¹Õµ•É¥Œ¡Ñ…É•Ð‘M¤°(€€€€€€€‘˜€ô…Ì¹¹Õµ•É¥Œ¡Ñ…É•Ð‘‘™}M…ÑÐ¤°(€€€€€€€Á}É…Ü€ô…Ì¹¹Õµ•É¥Œ¡Ñ…É•Ð‘Á}M…ÑÐ¤°(€€€€€€€Ý…É¹¥¹Ì€ôÁ…ÍÑ”¡Õ¹¥ÅÕ”¡…ÁÑÕÉ•‘}Ý…É¹¥¹Ì¤°½±±…ÁÍ”€ô€ˆð€ˆ¤°(€€€€€€€•ÉÉ½È€ô€ˆˆ(€€€€€€¤(€€€ô°Ý…É¹¥¹œ€ô™Õ¹Ñ¥½¸¡Ü¤ì(€€€€€…ÁÑÕÉ•‘}Ý…É¹¥¹Ì€ðð´Œ¡…ÁÑÕÉ•‘}Ý…É¹¥¹Ì°½¹‘¥Ñ¥½¹5•ÍÍ…”¡Ü¤¤(€€€€€¥¹Ù½­•I•ÍÑ…ÉÐ ‰µÕ™™±•]…É¹¥¹œˆ¤(€€€ô¤°(€€€•ÉÉ½È€ô™Õ¹Ñ¥½¸¡”¤ì(€€€€€±¥ÍÐ (€€€€€€€½¬€ô1M°(€€€€€€€‰•Ñ„€ô9}É•…±|°(€€€€€€€Í”€ô9}É•…±|°(€€€€€€€‘˜€ô9}É•…±|°(€€€€€€€Á}É…Ü€ô9}É•…±|°(€€€€€€€Ý…É¹¥¹Ì€ôÁ…ÍÑ”¡Õ¹¥ÅÕ”¡…ÁÑÕÉ•‘}Ý…É¹¥¹Ì¤°½±±…ÁÍ”€ô€ˆð€ˆ¤°(€€€€€€€•ÉÉ½È€ô½¹‘¥Ñ¥½¹5•ÍÍ…”¡”¤(€€€€€€¤(€€€ô(€€¤)ô()±½Í½}ÍÁ•Ì€ð´±¥ÍÐ (€±¥ÍÐ (€€€…¹…±åÍ¥Í}¥€ô€‰%I}‰±¥¹‘¥¹}É•Á½ÉÑ•ˆ°(€€€ÍÑÉ…ÑÕ´€ô€‰%¹Ñ•ÈµÉ•…‘•Èˆ°(€€€µ½‘•É…Ñ½È€ô€‰	±¥¹‘¥¹œÉ•Á½ÉÑ•ˆ°(€€€½¹ÑÉ…ÍÐ€ô€‰I•Á½ÉÑ•ÙÌ¹½ÐÉ•Á½ÉÑ•ˆ°(€€€‘…Ñ„€ô¥É}‰±¥¹‘¥¹œ°(€€€™½ÉµÕ±„€ôå¥}èø‰±¥¹‘¥¹}É•Á½ÉÑ•‘}‰¥¹…Éä°(€€€Ñ…É•Ñ}Ñ•É´€ô€‰‰±¥¹‘¥¹}É•Á½ÉÑ•‘}‰¥¹…Éäˆ°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì€ô9}¡…É…Ñ•É|(€€¤°(€±¥ÍÐ (€€€…¹…±åÍ¥Í}¥€ô€‰%I}É•¥½¹}ÕÉ½Á•}ÙÍ}9½ÉÑ¡}µ•É¥„ˆ°(€€€ÍÑÉ…ÑÕ´€ô€‰%¹Ñ•ÈµÉ•…‘•Èˆ°(€€€µ½‘•É…Ñ½È€ô€‰I•¥½¸ˆ°(€€€½¹ÑÉ…ÍÐ€ô€‰ÕÉ½Á”ÙÌ9½ÉÑ µ•É¥„ˆ°(€€€‘…Ñ„€ô¥É}É•¥½¸°(€€€™½ÉµÕ±„€ôå¥}èøÉ•¥½¹}™…Ñ½È°(€€€Ñ…É•Ñ}Ñ•É´€ô€‰É•¥½¹}™…Ñ½ÉÕÉ½Á”ˆ°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì€ô9}¡…É…Ñ•É|(€€¤°(€±¥ÍÐ (€€€…¹…±åÍ¥Í}¥€ô€‰%5}…Ñ•½É¥•Í|Ñ|Õ}ÙÍ|É|Ìˆ°(€€€ÍÑÉ…ÑÕ´€ô€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°(€€€µ½‘•É…Ñ½È€ô€‰9¼¸½˜…Ñ•½É¥•Ìˆ°(€€€½¹ÑÉ…ÍÐ€ô€ˆÐ´ÔÙÌ€È´Ìˆ°(€€€‘…Ñ„€ô¥µ}¹…Ð°(€€€™½ÉµÕ±„€ôå¥}èø¥Í}™Õ±°°(€€€Ñ…É•Ñ}Ñ•É´€ô€‰¥Í}™Õ±°ˆ°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì€ô9}¡…É…Ñ•É|(€€¤°(€±¥ÍÐ (€€€…¹…±åÍ¥Í}¥€ô€‰%5}µ½‘…±¥Ñå}Á…¥É}Q}UM}ÙÍ}Q}5I$ˆ°(€€€ÍÑÉ…ÑÕ´€ô€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°(€€€µ½‘•É…Ñ½È€ô€‰5½‘…±¥ÑäÁ…¥Èˆ°(€€€½¹ÑÉ…ÍÐ€ô€‰PµULµ‰…Í•ÙÌPµ5I$ˆ°(€€€‘…Ñ„€ô¥µ}Á…¥È°(€€€™½ÉµÕ±„€ôå¥}èøÁ…¥É}™…Ñ½È°(€€€Ñ…É•Ñ}Ñ•É´€ô€‰Á…¥É}™…Ñ½ÉQ}UM}ULˆ°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì€ô9}¡…É…Ñ•É|(€€¤°(€±¥ÍÐ (€€€…¹…±åÍ¥Í}¥€ô€‰%5}±½}Í…µÁ±•}Í¥é”ˆ°(€€€ÍÑÉ…ÑÕ´€ô€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°(€€€µ½‘•É…Ñ½È€ô€‰±½œ¡Í…µÁ±”Í¥é”¤ˆ°(€€€½¹ÑÉ…ÍÐ€ô€‰±½}¹Œˆ°(€€€‘…Ñ„€ô¥´€”ø”™¥±Ñ•È …¥Ì¹¹„¡±½}¹Œ¤¤°(€€€™½ÉµÕ±„€ôå¥}èø±½}¹Œ°(€€€Ñ…É•Ñ}Ñ•É´€ô€‰±½}¹Œˆ°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì€ô9}¡…É…Ñ•É|(€€¤°(€±¥ÍÐ (€€€…¹…±åÍ¥Í}¥€ô€‰%5}ÁÕ‰±¥…Ñ¥½¹}ÑåÁ•}…‰ÍÑÉ…Ñ}ÙÍ}™Õ±±Ñ•áÐˆ°(€€€ÍÑÉ…ÑÕ´€ô€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°(€€€µ½‘•É…Ñ½È€ô€‰AÕ‰±¥…Ñ¥½¸ÑåÁ”ˆ°(€€€½¹ÑÉ…ÍÐ€ô€‰‰ÍÑÉ…ÐÙÌ™Õ±°Ñ•áÐˆ°(€€€‘…Ñ„€ô¥µ}ÁÕ‰±¥…Ñ¥½¸°(€€€™½ÉµÕ±„€ôå¥}èø¥Í}…‰ÍÑÉ…Ð°(€€€Ñ…É•Ñ}Ñ•É´€ô€‰¥Í}…‰ÍÑÉ…Ðˆ°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì€ôÁ…ÍÑ” (€€€€€Í½ÉÐ¡Õ¹¥ÅÕ”¡¥µ}ÁÕ‰±¥…Ñ¥½¸‘ÍÑÕ‘å}¥‘m¥µ}ÁÕ‰±¥…Ñ¥½¸‘¥Í}…‰ÍÑÉ…Ð€ôô€Åt¤¤°(€€€€€½±±…ÁÍ”€ô€ˆìˆ(€€€€¤(€€¤°(€±¥ÍÐ (€€€…¹…±åÍ¥Í}¥€ô€‰%5}É•…‘•ÉÍ|ÍÁ±ÕÍ}ÙÍ|Äˆ°(€€€ÍÑÉ…ÑÕ´€ô€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°(€€€µ½‘•É…Ñ½È€ô€‰9¼¸½˜É•…‘•ÉÌˆ°(€€€½¹ÑÉ…ÍÐ€ô€ˆÌ½Èµ½É”ÙÌ€Äˆ°(€€€‘…Ñ„€ô¥µ}É•…‘•ÉÌ°(€€€™½ÉµÕ±„€ôå¥}èøÉ•…‘•É}½Õ¹Ñ}™…Ñ½È°(€€€Ñ…É•Ñ}Ñ•É´€ô€‰É•…‘•É}½Õ¹Ñ}™…Ñ½ÈÍÁ±ÕÌˆ°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì€ô9}¡…É…Ñ•É|(€€¤°(€±¥ÍÐ (€€€…¹…±åÍ¥Í}¥€ô€‰%5}É•ÁÉ•Í•¹Ñ…Ñ¥Ù•}ÍÁ•ÑÉÕ´ˆ°(€€€ÍÑÉ…ÑÕ´€ô€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ°(€€€µ½‘•É…Ñ½È€ô€‰I•ÁÉ•Í•¹Ñ…Ñ¥Ù”ÍÁ•ÑÉÕ´ˆ°(€€€½¹ÑÉ…ÍÐ€ô€‰e•ÌÙÌ¹¼½Õ¹±•…Èˆ°(€€€‘…Ñ„€ô¥µ}É•À°(€€€™½ÉµÕ±„€ôå¥}èøÉ•Á}å•Ì°(€€€Ñ…É•Ñ}Ñ•É´€ô€‰É•Á}å•Ìˆ°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì€ô9}¡…É…Ñ•É|(€€¤(¤()…¹½¹¥…±}¹½µ¥¹…°€ð´µ•Ñ…}É•É•ÍÍ¥½¹}É•ÍÕ±ÑÌ€”ø”(€™¥±Ñ•È¡Á}É…Ü€ð€À¸ÀÔ¤€”ø”(€Í•±•Ð¡ÍÑÉ…ÑÕ´°µ½‘•É…Ñ½È°½¹ÑÉ…ÍÐ°‰•Ñ„°Í”°‘˜°Á}É…Ü°‘™}¥¹Ñ•ÉÁÉ•Ñ…Ñ¥½¸¤()¥˜€¡¹É½Ü¡…¹½¹¥…±}¹½µ¥¹…°¤€„ô€à¤ì(€ÍÑ½À¡ÍÁÉ¥¹Ñ˜ (€€€€‰áÁ•Ñ••¥¡Ðµ½‘•É…Ñ½È½¹ÑÉ…ÍÑÌÝ¥Ñ ¹½µ¥¹…°ÀðÀ¸ÀÔ°™½Õ¹€•¸ˆ°(€€€¹É½Ü¡…¹½¹¥…±}¹½µ¥¹…°¤(€€¤¤)ô()ÍÁ•}­•åÌ€ð´‰¥¹‘}É½ÝÌ¡±…ÁÁ±ä¡±½Í½}ÍÁ•Ì°™Õ¹Ñ¥½¸¡à¤ì(€Ñ¥‰‰±” (€€€ÍÑÉ…ÑÕ´€ôà‘ÍÑÉ…ÑÕ´°(€€€µ½‘•É…Ñ½È€ôà‘µ½‘•É…Ñ½È°(€€€½¹ÑÉ…ÍÐ€ôà‘½¹ÑÉ…ÍÐ(€€¤)ô¤¤)­•å}½±Õµ¹Ì€ð´Œ ‰ÍÑÉ…ÑÕ´ˆ°€‰µ½‘•É…Ñ½Èˆ°€‰½¹ÑÉ…ÍÐˆ¤)¥˜€ (€¹É½Ü¡…¹Ñ¥}©½¥¸¡…¹½¹¥…±}¹½µ¥¹…°°ÍÁ•}­•åÌ°‰ä€ô­•å}½±Õµ¹Ì¤¤€ø€Àñð(€€€¹É½Ü¡…¹Ñ¥}©½¥¸¡ÍÁ•}­•åÌ°…¹½¹¥…±}¹½µ¥¹…°°‰ä€ô­•å}½±Õµ¹Ì¤¤€ø€À(¤ì(€ÍÑ½À ‰1=M<ÍÁ•¥™¥…Ñ¥½¹Ì‘¼¹½Ðµ…Ñ Ñ¡”¹½µ¥¹…°µÀµ•Ñ„µÉ•É•ÍÍ¥½¸É½ÝÌ¸ˆ¤)ô()±½Í½}Ù…±¥‘…Ñ¥½¹}Á…ÉÑÌ€ð´±¥ÍÐ ¤)±½Í½}¥Ñ•É…Ñ¥½¹}Á…ÉÑÌ€ð´±¥ÍÐ ¤()™½È€¡ÍÁ•Œ¥¸±½Í½}ÍÁ•Ì¤ì(€™Õ±±}™¥Ð€ð´™¥Ñ}±½Í½}Ñ…É•Ð¡ÍÁ•Œ‘‘…Ñ„°ÍÁ•Œ‘™½ÉµÕ±„°ÍÁ•Œ‘Ñ…É•Ñ}Ñ•É´¤(€¥˜€ …¥ÍQIU¡™Õ±±}™¥Ð‘½¬¤¤ì(€€€ÍÑ½À¡ÍÁÉ¥¹Ñ˜ (€€€€€€‰Õ±°µ½‘•°™…¥±•™½È€•Ìè€•Ìˆ°(€€€€€ÍÁ•Œ‘…¹…±åÍ¥Í}¥°(€€€€€™Õ±±}™¥Ð‘•ÉÉ½È(€€€€¤¤(€ô((€…¹½¹¥…±}É½Ü€ð´…¹½¹¥…±}¹½µ¥¹…°€”ø”(€€€™¥±Ñ•È (€€€€€ÍÑÉ…ÑÕ´€ôôÍÁ•Œ‘ÍÑÉ…ÑÕ´°(€€€€€µ½‘•É…Ñ½È€ôôÍÁ•Œ‘µ½‘•É…Ñ½È°(€€€€€½¹ÑÉ…ÍÐ€ôôÍÁ•Œ‘½¹ÑÉ…ÍÐ(€€€€¤(€¥˜€¡¹É½Ü¡…¹½¹¥…±}É½Ü¤€„ô€Ä¤ì(€€€ÍÑ½À¡ÍÁÉ¥¹Ñ˜ ‰…¹½¹¥…°É½Ü±½½­ÕÀ™…¥±•™½È€•Ì¸ˆ°ÍÁ•Œ‘…¹…±åÍ¥Í}¥¤¤(€ô((€±½Í½}Ù…±¥‘…Ñ¥½¹}Á…ÉÑÍmm±•¹Ñ ¡±½Í½}Ù…±¥‘…Ñ¥½¹}Á…ÉÑÌ¤€¬€Åut€ð´Ñ¥‰‰±” (€€€…¹…±åÍ¥Í}¥€ôÍÁ•Œ‘…¹…±åÍ¥Í}¥°(€€€ÍÑÉ…ÑÕ´€ôÍÁ•Œ‘ÍÑÉ…ÑÕ´°(€€€µ½‘•É…Ñ½È€ôÍÁ•Œ‘µ½‘•É…Ñ½È°(€€€½¹ÑÉ…ÍÐ€ôÍÁ•Œ‘½¹ÑÉ…ÍÐ°(€€€Ñ…É•Ñ}Ñ•É´€ôÍÁ•Œ‘Ñ…É•Ñ}Ñ•É´°(€€€•™™•ÑÌ€ô¹É½Ü¡ÍÁ•Œ‘‘…Ñ„¤°(€€€ÍÑÕ‘¥•Ì€ô¹}‘¥ÍÑ¥¹Ð¡ÍÁ•Œ‘‘…Ñ„‘ÍÑÕ‘å}¥¤°(€€€±ÕÍÑ•ÉÌ€ô¹}‘¥ÍÑ¥¹Ð¡ÍÁ•Œ‘‘…Ñ„‘‘•Á}¥¤°(€€€‰•Ñ…}É•½µÁÕÑ•€ô™Õ±±}™¥Ð‘‰•Ñ„°(€€€‰•Ñ…}…¹½¹¥…°€ô…¹½¹¥…±}É½Ü‘‰•Ñ„°(€€€‰•Ñ…}…‰Í}‘¥™˜€ô…‰Ì¡™Õ±±}™¥Ð‘‰•Ñ„€´…¹½¹¥…±}É½Ü‘‰•Ñ„¤°(€€€Í•}É•½µÁÕÑ•€ô™Õ±±}™¥Ð‘Í”°(€€€Í•}…¹½¹¥…°€ô…¹½¹¥…±}É½Ü‘Í”°(€€€Í•}…‰Í}‘¥™˜€ô…‰Ì¡™Õ±±}™¥Ð‘Í”€´…¹½¹¥…±}É½Ü‘Í”¤°(€€€‘™}É•½µÁÕÑ•€ô™Õ±±}™¥Ð‘‘˜°(€€€‘™}…¹½¹¥…°€ô…¹½¹¥…±}É½Ü‘‘˜°(€€€‘™}…‰Í}‘¥™˜€ô…‰Ì¡™Õ±±}™¥Ð‘‘˜€´…¹½¹¥…±}É½Ü‘‘˜¤°(€€€Á}É•½µÁÕÑ•€ô™Õ±±}™¥Ð‘Á}É…Ü°(€€€Á}…¹½¹¥…°€ô…¹½¹¥…±}É½Ü‘Á}É…Ü°(€€€Á}…‰Í}‘¥™˜€ô…‰Ì¡™Õ±±}™¥Ð‘Á}É…Ü€´…¹½¹¥…±}É½Ü‘Á}É…Ü¤°(€€€…¹½¹¥…±}¥¹Ñ•ÉÁÉ•Ñ…Ñ¥½¸€ô…¹½¹¥…±}É½Ü‘‘™}¥¹Ñ•ÉÁÉ•Ñ…Ñ¥½¸°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì€ôÍÁ•Œ‘Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì°(€€€™Õ±±}µ½‘•±}Ý…É¹¥¹Ì€ô™Õ±±}™¥Ð‘Ý…É¹¥¹Ì(€€¤((€™½È€¡½µ¥ÑÑ•‘}ÍÑÕ‘ä¥¸Í½ÉÐ¡Õ¹¥ÅÕ”¡ÍÁ•Œ‘‘…Ñ„‘ÍÑÕ‘å}¥¤¤¤ì(€€€‘•±•Ñ¥½¹}‘…Ñ„€ð´ÍÁ•Œ‘‘…Ñ„€”ø”™¥±Ñ•È¡ÍÑÕ‘å}¥€„ô½µ¥ÑÑ•‘}ÍÑÕ‘ä¤(€€€‘•±•Ñ¥½¹}™¥Ð€ð´™¥Ñ}±½Í½}Ñ…É•Ð (€€€€€‘•±•Ñ¥½¹}‘…Ñ„°(€€€€€ÍÁ•Œ‘™½ÉµÕ±„°(€€€€€ÍÁ•Œ‘Ñ…É•Ñ}Ñ•É´(€€€€¤(€€€±½Í½}¥Ñ•É…Ñ¥½¹}Á…ÉÑÍmm±•¹Ñ ¡±½Í½}¥Ñ•É…Ñ¥½¹}Á…ÉÑÌ¤€¬€Åut€ð´Ñ¥‰‰±” (€€€€€…¹…±åÍ¥Í}¥€ôÍÁ•Œ‘…¹…±åÍ¥Í}¥°(€€€€€ÍÑÉ…ÑÕ´€ôÍÁ•Œ‘ÍÑÉ…ÑÕ´°(€€€€€µ½‘•É…Ñ½È€ôÍÁ•Œ‘µ½‘•É…Ñ½È°(€€€€€½¹ÑÉ…ÍÐ€ôÍÁ•Œ‘½¹ÑÉ…ÍÐ°(€€€€€½µ¥ÑÑ•‘}ÍÑÕ‘ä€ô½µ¥ÑÑ•‘}ÍÑÕ‘ä°(€€€€€É•µ½Ù•‘}•™™•ÑÌ€ôÍÕ´¡ÍÁ•Œ‘‘…Ñ„‘ÍÑÕ‘å}¥€ôô½µ¥ÑÑ•‘}ÍÑÕ‘ä¤°(€€€€€É•µ…¥¹¥¹}•™™•ÑÌ€ô¹É½Ü¡‘•±•Ñ¥½¹}‘…Ñ„¤°(€€€€€É•µ…¥¹¥¹}ÍÑÕ‘¥•Ì€ô¹}‘¥ÍÑ¥¹Ð¡‘•±•Ñ¥½¹}‘…Ñ„‘ÍÑÕ‘å}¥¤°(€€€€€É•µ…¥¹¥¹}±ÕÍÑ•ÉÌ€ô¹}‘¥ÍÑ¥¹Ð¡‘•±•Ñ¥½¹}‘…Ñ„‘‘•Á}¥¤°(€€€€€ÍÕ•ÍÍ™Õ°€ô‘•±•Ñ¥½¹}™¥Ð‘½¬°(€€€€€‰•Ñ„€ô‘•±•Ñ¥½¹}™¥Ð‘‰•Ñ„°(€€€€€Í”€ô‘•±•Ñ¥½¹}™¥Ð‘Í”°(€€€€€‘˜€ô‘•±•Ñ¥½¹}™¥Ð‘‘˜°(€€€€€Á}É…Ü€ô‘•±•Ñ¥½¹}™¥Ð‘Á}É…Ü°(€€€€€¹½µ¥¹…±}Á}±Ñ|Á|ÀÔ€ô¥™•±Í” (€€€€€€€‘•±•Ñ¥½¹}™¥Ð‘½¬°(€€€€€€€‘•±•Ñ¥½¹}™¥Ð‘Á}É…Ü€ð€À¸ÀÔ°(€€€€€€€9(€€€€€€¤°(€€€€€‘¥É•Ñ¥½¹}½¹Í¥ÍÑ•¹Ð€ô¥™•±Í” (€€€€€€€‘•±•Ñ¥½¹}™¥Ð‘½¬°(€€€€€€€Í¥¸¡‘•±•Ñ¥½¹}™¥Ð‘‰•Ñ„¤€ôôÍ¥¸¡™Õ±±}™¥Ð‘‰•Ñ„¤°(€€€€€€€9(€€€€€€¤°(€€€€€Ý…É¹¥¹Ì€ô‘•±•Ñ¥½¹}™¥Ð‘Ý…É¹¥¹Ì°(€€€€€™…¥±ÕÉ•}É•…Í½¸€ô‘•±•Ñ¥½¹}™¥Ð‘•ÉÉ½È(€€€€¤(€ô)ô()±½Í½}Ù…±¥‘…Ñ¥½¸€ð´‰¥¹‘}É½ÝÌ¡±½Í½}Ù…±¥‘…Ñ¥½¹}Á…ÉÑÌ¤)±½Í½}¥Ñ•É…Ñ¥½¹Ì€ð´‰¥¹‘}É½ÝÌ¡±½Í½}¥Ñ•É…Ñ¥½¹}Á…ÉÑÌ¤()Ù…±¥‘…Ñ¥½¹}Ñ½±•É…¹”€ð´€Å”´ÄÀ)¥˜€¡…¹ä (€±½Í½}Ù…±¥‘…Ñ¥½¸‘‰•Ñ…}…‰Í}‘¥™˜€øÙ…±¥‘…Ñ¥½¹}Ñ½±•É…¹”ð(€€€±½Í½}Ù…±¥‘…Ñ¥½¸‘Í•}…‰Í}‘¥™˜€øÙ…±¥‘…Ñ¥½¹}Ñ½±•É…¹”ð(€€€±½Í½}Ù…±¥‘…Ñ¥½¸‘‘™}…‰Í}‘¥™˜€øÙ…±¥‘…Ñ¥½¹}Ñ½±•É…¹”ð(€€€±½Í½}Ù…±¥‘…Ñ¥½¸‘Á}…‰Í}‘¥™˜€øÙ…±¥‘…Ñ¥½¹}Ñ½±•É…¹”(¤¤ì(€ÍÑ½À ‰Ð±•…ÍÐ½¹”1=M<™Õ±°µ½‘•°‘¥¹½Ðµ…Ñ Ñ¡”…¹½¹¥…°É•ÍÕ±Ð¸ˆ¤)ô()±½Í½}ÍÕµµ…Éä€ð´±½Í½}¥Ñ•É…Ñ¥½¹Ì€”ø”(€±•™Ñ}©½¥¸ (€€€±½Í½}Ù…±¥‘…Ñ¥½¸€”ø”(€€€€€Í•±•Ð (€€€€€€€…¹…±åÍ¥Í}¥°(€€€€€€€™Õ±±}‰•Ñ„€ô‰•Ñ…}É•½µÁÕÑ•°(€€€€€€€™Õ±±}Á}É…Ü€ôÁ}É•½µÁÕÑ•°(€€€€€€€…¹½¹¥…±}¥¹Ñ•ÉÁÉ•Ñ…Ñ¥½¸°(€€€€€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì(€€€€€€¤°(€€€‰ä€ô€‰…¹…±åÍ¥Í}¥ˆ(€€¤€”ø”(€É½ÕÁ}‰ä (€€€…¹…±åÍ¥Í}¥°ÍÑÉ…ÑÕ´°µ½‘•É…Ñ½È°½¹ÑÉ…ÍÐ°(€€€™Õ±±}‰•Ñ„°™Õ±±}Á}É…Ü°…¹½¹¥…±}¥¹Ñ•ÉÁÉ•Ñ…Ñ¥½¸°(€€€Ñ…É•Ñ}±•Ù•±}ÍÑÕ‘¥•Ì(€€¤€”ø”(€ÍÕµµ…É¥Í” (€€€Á±…¹¹•‘}¥Ñ•É…Ñ¥½¹Ì€ô¸ ¤°(€€€ÍÕ•ÍÍ™Õ±}¥Ñ•É…Ñ¥½¹Ì€ôÍÕ´¡ÍÕ•ÍÍ™Õ°¤°(€€€¹½¹}•ÍÑ¥µ…‰±•}¥Ñ•É…Ñ¥½¹Ì€ôÍÕ´ …ÍÕ•ÍÍ™Õ°¤°(€€€‘¥É•Ñ¥½¹}½¹Í¥ÍÑ•¹Ñ}¸€ôÍÕ´¡‘¥É•Ñ¥½¹}½¹Í¥ÍÑ•¹Ð°¹„¹É´€ôQIU¤°(€€€‘¥É•Ñ¥½¹}½¹Í¥ÍÑ•¹Ñ}‘•¹½µ¥¹…Ñ½È€ôÍÕ´¡ÍÕ•ÍÍ™Õ°¤°(€€€¹½µ¥¹…±}Á}±Ñ|Á|ÀÕ}¸€ôÍÕ´¡¹½µ¥¹…±}Á}±Ñ|Á|ÀÔ°¹„¹É´€ôQIU¤°(€€€¹½µ¥¹…±}Á}±Ñ|Á|ÀÕ}‘•¹½µ¥¹…Ñ½È€ôÍÕ´¡ÍÕ•ÍÍ™Õ°¤°(€€€‰•Ñ…}µ¥¸€ô¥™•±Í”¡…¹ä¡ÍÕ•ÍÍ™Õ°¤°µ¥¸¡‰•Ñ…mÍÕ•ÍÍ™Õ±t¤°9}É•…±|¤°(€€€‰•Ñ…}µ…à€ô¥™•±Í”¡…¹ä¡ÍÕ•ÍÍ™Õ°¤°µ…à¡‰•Ñ…mÍÕ•ÍÍ™Õ±t¤°9}É•…±|¤°(€€€Á}É…Ý}µ¥¸€ô¥™•±Í”¡…¹ä¡ÍÕ•ÍÍ™Õ°¤°µ¥¸¡Á}É…ÝmÍÕ•ÍÍ™Õ±t¤°9}É•…±|¤°(€€€Á}É…Ý}µ…à€ô¥™•±Í”¡…¹ä¡ÍÕ•ÍÍ™Õ°¤°µ…à¡Á}É…ÝmÍÕ•ÍÍ™Õ±t¤°9}É•…±|¤°(€€€¹½¹}•ÍÑ¥µ…‰±•}ÍÑÕ‘¥•Ì€ôÁ…ÍÑ”¡½µ¥ÑÑ•‘}ÍÑÕ‘ål…ÍÕ•ÍÍ™Õ±t°½±±…ÁÍ”€ô€ˆìˆ¤°(€€€™…¥±ÕÉ•}É•…Í½¹Ì€ôÁ…ÍÑ”¡Õ¹¥ÅÕ”¡™…¥±ÕÉ•}É•…Í½¹l…ÍÕ•ÍÍ™Õ±t¤°½±±…ÁÍ”€ô€ˆð€ˆ¤°(€€€Ý…É¹¥¹}¥Ñ•É…Ñ¥½¹Ì€ôÍÕ´¡Ý…É¹¥¹Ì€„ô€ˆˆ¤°(€€€€¹É½ÕÁÌ€ô€‰‘É½Àˆ(€€¤()ÝÉ¥Ñ•}ÍØ (€±½Í½}¥Ñ•É…Ñ¥½¹Ì°(€™¥±”¹Á…Ñ ¡½ÕÑ}‘¥È°€‰±½Í½}¹½µ¥¹…±}µ½‘•É…Ñ½É}¥Ñ•É…Ñ¥½¹Ì¹ÍØˆ¤(¤)ÝÉ¥Ñ•}ÍØ (€±½Í½}ÍÕµµ…Éä°(€™¥±”¹Á…Ñ ¡½ÕÑ}‘¥È°€‰±½Í½}¹½µ¥¹…±}µ½‘•É…Ñ½É}ÍÕµµ…Éä¹ÍØˆ¤(¤)ÝÉ¥Ñ•}ÍØ (€±½Í½}Ù…±¥‘…Ñ¥½¸°(€™¥±”¹Á…Ñ ¡½ÕÑ}‘¥È°€‰±½Í½}¹½µ¥¹…±}µ½‘•É…Ñ½É}Ù…±¥‘…Ñ¥½¸¹ÍØˆ¤(¤((ŒŒ½ÕÍ•±•…Ù”µ½¹”µ…‰ÍÑÉ…ÐµÍÑÕ‘äµ½ÕÐ¡•­Ì™½ÈÑ¡”¥¹Ñ•ÈµÉ•…‘•È(ŒŒÁÕ‰±¥…Ñ¥½¸µÑåÁ”½¹ÑÉ…ÍÐ¸Q¡•Í”…É”Í•Á…É…Ñ”™É½´Ñ¡”¹½µ¥¹…°µÀÍ•Ð…‰½Ù”(ŒŒ‰•…ÕÍ”Ñ¡”™Õ±°¥¹Ñ•ÈµÉ•…‘•ÈÁÕ‰±¥…Ñ¥½¸µÑåÁ”½¹ÑÉ…ÍÐ¡…ÀôÀ¸ÈÜÈ¸()ÁÕ‰±¥…Ñ¥½¹}¡•­}Á…ÉÑÌ€ð´±¥ÍÐ ¤)ÁÕ‰±¥…Ñ¥½¹}¡•­}ÍÁ•Ì€ð´Œ ‰Õ±°µ½‘•°ˆ°Í½ÉÐ¡Õ¹¥ÅÕ” (€¥É}ÁÕ‰±¥…Ñ¥½¸‘ÍÑÕ‘å}¥‘m¥É}ÁÕ‰±¥…Ñ¥½¸‘¥Í}…‰ÍÑÉ…Ð€ôô€Åt(¤¤¤()™½È€¡¡•­}±…‰•°¥¸ÁÕ‰±¥…Ñ¥½¹}¡•­}ÍÁ•Ì¤ì(€½µ¥ÑÑ•‘}ÍÑÕ‘ä€ð´¥˜€¡¡•­}±…‰•°€ôô€‰Õ±°µ½‘•°ˆ¤9}¡…É…Ñ•É|•±Í”¡•­}±…‰•°(€¡•­}‘…Ñ„€ð´¥˜€¡¥Ì¹¹„¡½µ¥ÑÑ•‘}ÍÑÕ‘ä¤¤ì(€€€¥É}ÁÕ‰±¥…Ñ¥½¸(€ô•±Í”ì(€€€¥É}ÁÕ‰±¥…Ñ¥½¸€”ø”™¥±Ñ•È¡ÍÑÕ‘å}¥€„ô½µ¥ÑÑ•‘}ÍÑÕ‘ä¤(€ô(€¡•­}™¥Ð€ð´™¥Ñ}±½Í½}Ñ…É•Ð (€€€¡•­}‘…Ñ„°(€€€å¥}èø¥Í}…‰ÍÑÉ…Ð°(€€€€‰¥Í}…‰ÍÑÉ…Ðˆ(€€¤(€…‰ÍÑÉ…Ñ}‘…Ñ„€ð´¡•­}‘…Ñ„€”ø”™¥±Ñ•È¡¥Í}…‰ÍÑÉ…Ð€ôô€Ä¤(€…‰ÍÑÉ…Ñ}ÍÑÕ‘¥•Ì€ð´¹}‘¥ÍÑ¥¹Ð¡…‰ÍÑÉ…Ñ}‘…Ñ„‘ÍÑÕ‘å}¥¤(€…‰ÍÑÉ…Ñ}±ÕÍÑ•ÉÌ€ð´¹}‘¥ÍÑ¥¹Ð¡…‰ÍÑÉ…Ñ}‘…Ñ„‘‘•Á}¥¤(€ÁÕ‰±¥…Ñ¥½¹}¡•­}Á…ÉÑÍmm±•¹Ñ ¡ÁÕ‰±¥…Ñ¥½¹}¡•­}Á…ÉÑÌ¤€¬€Åut€ð´Ñ¥‰‰±” (€€€…¹…±åÍ¥Ì€ô¡•­}±…‰•°°(€€€½µ¥ÑÑ•‘}ÍÑÕ‘ä€ô½µ¥ÑÑ•‘}ÍÑÕ‘ä°(€€€•™™•ÑÌ€ô¹É½Ü¡¡•­}‘…Ñ„¤°(€€€ÍÑÕ‘¥•Ì€ô¹}‘¥ÍÑ¥¹Ð¡¡•­}‘…Ñ„‘ÍÑÕ‘å}¥¤°(€€€±ÕÍÑ•ÉÌ€ô¹}‘¥ÍÑ¥¹Ð¡¡•­}‘…Ñ„‘‘•Á}¥¤°(€€€…‰ÍÑÉ…Ñ}•™™•ÑÌ€ô¹É½Ü¡…‰ÍÑÉ…Ñ}‘…Ñ„¤°(€€€…‰ÍÑÉ…Ñ}ÍÑÕ‘¥•Ì€ô…‰ÍÑÉ…Ñ}ÍÑÕ‘¥•Ì°(€€€…‰ÍÑÉ…Ñ}±ÕÍÑ•ÉÌ€ô…‰ÍÑÉ…Ñ}±ÕÍÑ•ÉÌ°(€€€ÍÕ•ÍÍ™Õ°€ô¡•­}™¥Ð‘½¬°(€€€‰•Ñ„€ô¡•­}™¥Ð‘‰•Ñ„°(€€€Í”€ô¡•­}™¥Ð‘Í”°(€€€‘˜€ô¡•­}™¥Ð‘‘˜°(€€€Á}É…Ü€ô¡•­}™¥Ð‘Á}É…Ü°(€€€¥¹Ñ•ÉÁÉ•Ñ…Ñ¥½¸€ô…Í•}Ý¡•¸ (€€€€€€…¡•­}™¥Ð‘½¬ø€‰9½¸µ•ÍÑ¥µ…‰±”ˆ°(€€€€€…‰ÍÑÉ…Ñ}ÍÑÕ‘¥•Ì€ð€Èñð…‰ÍÑÉ…Ñ}±ÕÍÑ•ÉÌ€ð€Èø(€€€€€€€€‰9½¸µ¥¹Ñ•ÉÁÉ•Ñ…‰±”€¡Í¥¹±•Ñ½¸…‰ÍÑÉ…Ð±•Ù•°¤ˆ°(€€€€€¡•­}™¥Ð‘‘˜€ð€Ðø€‰9½¸µ¥¹Ñ•ÉÁÉ•Ñ…‰±”€¡‘˜ðÐ¤ˆ°(€€€€€¡•­}™¥Ð‘‘˜€ð€ÄÀø€‰É…¥±”€¡‘˜€Ð´ðÄÀ¤ˆ°(€€€€€QIUø€‰áÁ±½É…Ñ½Éäˆ(€€€€¤°(€€€Ý…É¹¥¹Ì€ô¡•­}™¥Ð‘Ý…É¹¥¹Ì°(€€€™…¥±ÕÉ•}É•…Í½¸€ô¡•­}™¥Ð‘•ÉÉ½È(€€¤)ô()¥¹Ñ•É}É•…‘•É}ÁÕ‰±¥…Ñ¥½¹}¡•­Ì€ð´‰¥¹‘}É½ÝÌ¡ÁÕ‰±¥…Ñ¥½¹}¡•­}Á…ÉÑÌ¤)ÝÉ¥Ñ•}ÍØ (€¥¹Ñ•É}É•…‘•É}ÁÕ‰±¥…Ñ¥½¹}¡•­Ì°(€™¥±”¹Á…Ñ ¡½ÕÑ}‘¥È°€‰¥¹Ñ•É}É•…‘•É}ÁÕ‰±¥…Ñ¥½¹}ÑåÁ•}±•…Ù•}½¹•}…‰ÍÑÉ…Ð¹ÍØˆ¤(¤()…Ð ‰q¹9½µ¥¹…°µÀµ½‘•É…Ñ½È±•…Ù”µ½¹”µÍÑÕ‘äµ½ÕÐÍÕµµ…Éåq¸ˆ¤)ÁÉ¥¹Ð¡±½Í½}ÍÕµµ…Éä°¸€ô%¹˜°Ý¥‘Ñ €ô%¹˜¤)…Ð ‰q¹%¹Ñ•ÈµÉ•…‘•ÈÁÕ‰±¥…Ñ¥½¸µÑåÁ”±•…Ù”µ½¹”µ…‰ÍÑÉ…Ð¡•­Íq¸ˆ¤)ÁÉ¥¹Ð¡¥¹Ñ•É}É•…‘•É}ÁÕ‰±¥…Ñ¥½¹}¡•­Ì°¸€ô%¹˜°Ý¥‘Ñ €ô%¹˜¤