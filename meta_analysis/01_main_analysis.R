## ============================================================
## Bosniak Classification Reproducibility SRMA
## Main Meta-Analysis — Revised Framework
##
## Primary: kappa-only, RVE+CR2
## Sensitivity: weighted-only, unweighted-only, 5-cat-only
## Exploratory: all metrics pooled
## ============================================================

library(metafor)
library(clubSandwich)
library(dplyr)
library(readr)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 0) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
setwd(script_dir)

# ============================================================
# 1. DATA LOADING & PREPARATION
# ============================================================

comp <- read_csv("comparisons.csv", show_col_types = FALSE)
stud <- read_csv("studies.csv", show_col_types = FALSE)
qarel <- read_csv("qarel_for_analysis.csv", show_col_types = FALSE)

# Publication type is maintained in studies.csv as the single source of truth.
standard_publication_types <- c("journal_article", "conference_abstract")
active_studies <- unique(comp$study_id[is.na(comp$exclude) | comp$exclude == ""])
active_publication_types <- stud %>%
  filter(study_id %in% active_studies) %>%
  pull(publication_type)
if (any(is.na(active_publication_types)) ||
    any(!active_publication_types %in% standard_publication_types)) {
  stop("Every retained study must have a standardized publication_type in studies.csv")
}
if (stud$publication_type[stud$study_id == "Chang_2015"] != "conference_abstract") {
  stop("Chang_2015 must be coded as conference_abstract in studies.csv")
}

# Join study-level variables
comp <- comp %>%
  left_join(stud %>% select(study_id, year, publication_type, study_design),
            by = "study_id") %>%
  left_join(qarel %>% select(study_id, overall_rob, qarel_yes_count),
            by = "study_id")

dat <- comp %>%
  filter(is.na(exclude) | exclude == "") %>%
  mutate(
    kappa = as.numeric(kappa),
    gwet_ac2 = as.numeric(gwet_ac2),
    icc = as.numeric(icc),
    ci_lo = as.numeric(ci_lower),
    ci_hi = as.numeric(ci_upper),
    se_raw = as.numeric(kappa_se),
    nc = as.numeric(n_cysts),
    n_readers = as.numeric(n_readers),
    n_categories = as.numeric(n_categories),
    stat = coalesce(kappa, gwet_ac2, icc),
    stat_type = case_when(
      !is.na(kappa) ~ "kappa",
      !is.na(gwet_ac2) ~ "gwet",
      !is.na(icc) ~ "icc",
      TRUE ~ NA_character_
    ),

    # ---- Moderator variable coding ----

    # Weight: classified (weighted/unweighted) for primary; 3-level for sensitivity
    wt_group = case_when(
      weight_scheme %in% c("linear", "quadratic", "custom",
                           "Cicchetti", "weighted_unspecified") ~ "weighted",
      weight_scheme == "unweighted" ~ "unweighted",
      TRUE ~ "other"
    ),
    wt_3level = case_when(
      weight_scheme %in% c("linear", "quadratic", "custom",
                           "Cicchetti", "weighted_unspecified") ~ "weighted",
      weight_scheme == "unweighted" ~ "unweighted",
      weight_scheme == "unspecified" ~ "unspecified",
      TRUE ~ NA_character_
    ),

    # n_categories: 2-3 vs 4-5
    ncat_group = case_when(
      n_categories <= 3 ~ "2-3",
      n_categories <= 5 ~ "4-5",
      n_categories > 5 ~ "6+",
      TRUE ~ NA_character_
    ),

    # Modality: CT vs MRI vs CEUS_US vs other
    mod_group = case_when(
      std_modality_1 == "CT" ~ "CT",
      std_modality_1 == "MRI" ~ "MRI",
      std_modality_1 %in% c("CEUS", "US", "SMI") ~ "CEUS_US",
      TRUE ~ "other"
    ),

    # Version: v2019 vs pre-2019
    ver_group = case_when(
      std_version_1 == "v2019" ~ "v2019",
      std_version_1 %in% c("original", "binary", "EFSUMB", "modified") ~ "pre2019",
      TRUE ~ NA_character_
    ),

    # Publication type
    pub_type = case_when(
      publication_type == "journal_article" ~ "fulltext",
      publication_type == "conference_abstract" ~ "abstract",
      TRUE ~ NA_character_
    ),

    # n_readers: 2 vs 3+
    nr_group = case_when(
      n_readers == 2 ~ "2",
      n_readers >= 3 ~ "3plus",
      TRUE ~ NA_character_
    ),

    # log(n_cysts) for continuous moderator
    log_nc = ifelse(!is.na(nc) & nc > 0, log(nc), NA_real_),

    # Reader experience (continuous, parsed mean years)
    exp_mean = as.numeric(ifelse(exp_mean == "", NA, exp_mean)),

    # Specialty: subspecialty vs general (binary)
    spec_group = case_when(
      specialty_group == "subspecialty" ~ "subspecialty",
      specialty_group %in% c("general", "mixed", "residents",
                             "non_radiology") ~ "non_subspecialty",
      TRUE ~ NA_character_
    ),

    # Blinding: blinded vs not (binary)
    is_blinded = case_when(
      blinding_reported == "blinded" ~ 1L,
      blinding_reported %in% c("unblinded", "unclear") ~ 0L,
      TRUE ~ NA_integer_
    ),

    # Reader structure: pooled vs non-pooled (binary)
    is_pooled = case_when(
      reader_structure == "pooled" ~ 1L,
      reader_structure %in% c("single_pair", "pairwise",
                              "single", "subgroup") ~ 0L,
      TRUE ~ NA_integer_
    ),

    # Study design: prospective vs retrospective
    is_prospective = case_when(
      study_design_cat == "prospective" ~ 1L,
      study_design_cat == "retrospective" ~ 0L,
      TRUE ~ NA_integer_
    ),

    # Region (3-level: North_America, Asia, Europe; Other/empty → NA)
    region_cat = case_when(
      region %in% c("North_America", "Asia", "Europe") ~ region,
      TRUE ~ NA_character_
    ),

    # Publication year (continuous, from join)
    pub_year = as.numeric(year),

    # QAREL quality: Low vs Moderate/High (binary)
    is_low_rob = case_when(
      overall_rob == "Low" ~ 1L,
      overall_rob %in% c("Moderate", "High") ~ 0L,
      TRUE ~ NA_integer_
    ),
    # QAREL yes count (continuous)
    qarel_score = as.numeric(qarel_yes_count),

    # Modality pair (for inter-modality)
    mod_pair = case_when(
      (std_modality_1 %in% c("CT","DECT") & std_modality_2 == "MRI") |
        (std_modality_1 == "MRI" & std_modality_2 %in% c("CT","DECT")) ~ "CT_MRI",
      (std_modality_1 %in% c("CT","CT_MRI") & std_modality_2 %in% c("CEUS","US")) |
        (std_modality_1 %in% c("CEUS","US") & std_modality_2 %in% c("CT","CT_MRI")) ~ "CT_CEUS",
      (std_modality_1 == "MRI" & std_modality_2 %in% c("CEUS","US")) |
        (std_modality_1 %in% c("CEUS","US") & std_modality_2 == "MRI") ~ "MRI_CEUS",
      TRUE ~ "other"
    )
  ) %>%
  filter(!is.na(stat))

# Fisher's z + variance
dat <- dat %>%
  mutate(
    stat_clamped = pmin(pmax(stat, -0.995), 0.995),
    yi = atanh(stat_clamped),
    vi_ci = ifelse(!is.na(ci_lo) & !is.na(ci_hi),
                   ((atanh(pmin(pmax(ci_hi, -0.995), 0.995)) -
                     atanh(pmin(pmax(ci_lo, -0.995), 0.995))) / (2 * qnorm(0.975)))^2,
                   NA_real_),
    vi_se = ifelse(!is.na(se_raw),
                   (se_raw / (1 - stat_clamped^2))^2, NA_real_),
    vi_n = ifelse(!is.na(nc) & nc > 3, 1 / (nc - 3), NA_real_),
    vi = coalesce(vi_ci, vi_se, vi_n)
  ) %>%
  filter(!is.na(vi) & vi > 0)

cat(sprintf("Total valid effects: %d\n\n", nrow(dat)))


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

run_rve_cr2 <- function(data, label = "") {
  if (nrow(data) < 3) {
    cat(sprintf("\n--- %s: SKIPPED (k=%d < 3) ---\n", label, nrow(data)))
    return(NULL)
  }
  n_clusters <- length(unique(data$dep_id))
  if (n_clusters < 3) {
    cat(sprintf("\n--- %s: SKIPPED (clusters=%d < 3) ---\n", label, n_clusters))
    return(NULL)
  }

  mod <- rma.mv(yi, vi,
                random = list(~ 1 | dep_id, ~ 1 | comparison_id),
                data = data, method = "REML")
  cr2 <- coef_test(mod, vcov = "CR2", cluster = data$dep_id)

  pooled_z <- mod$b[1]
  pooled_kappa <- tanh(pooled_z)
  cr2_se <- cr2$SE; cr2_df <- cr2$df_Satt
  t_crit <- qt(0.975, cr2_df)
  ci_kappa <- tanh(c(pooled_z - t_crit * cr2_se, pooled_z + t_crit * cr2_se))

  pi_z <- c(pooled_z - t_crit * sqrt(cr2_se^2 + mod$sigma2[1]),
            pooled_z + t_crit * sqrt(cr2_se^2 + mod$sigma2[1]))
  pi_kappa <- tanh(pi_z)

  k <- nrow(data)
  W <- 1 / data$vi
  typical_v <- (k - 1) * sum(W) / (sum(W)^2 - sum(W^2))
  I2_total <- 100 * sum(mod$sigma2) / (sum(mod$sigma2) + typical_v)
  I2_between <- 100 * mod$sigma2[1] / (sum(mod$sigma2) + typical_v)
  I2_within <- 100 * mod$sigma2[2] / (sum(mod$sigma2) + typical_v)

  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("Effects: %d, Clusters: %d, Studies: %d\n",
              k, n_clusters, length(unique(data$study_id))))
  cat(sprintf("Pooled: %.3f (95%% CI: %.3f, %.3f)\n",
              pooled_kappa, ci_kappa[1], ci_kappa[2]))
  cat(sprintf("95%% PI: %.3f to %.3f\n", pi_kappa[1], pi_kappa[2]))
  cat(sprintf("CR2: SE=%.4f, df=%.1f, p=%s\n",
              cr2_se, cr2_df, format.pval(cr2$p_Satt, digits = 3)))
  cat(sprintf("tau2: between=%.4f, within=%.4f\n", mod$sigma2[1], mod$sigma2[2]))
  cat(sprintf("I2: total=%.1f%%, between=%.1f%%, within=%.1f%%\n",
              I2_total, I2_between, I2_within))

  list(mod = mod, cr2 = cr2, pooled_kappa = pooled_kappa,
       ci_kappa = ci_kappa, pi_kappa = pi_kappa,
       I2_total = I2_total, I2_between = I2_between, I2_within = I2_within,
       k = k, n_clusters = n_clusters, cr2_df = cr2_df)
}

run_opes_re <- function(data, label = "") {
  if (nrow(data) < 3) {
    cat(sprintf("\n--- %s (OPES): SKIPPED (k=%d) ---\n", label, nrow(data)))
    return(NULL)
  }
  data <- data %>%
    mutate(
      opes_val_str = ifelse(
        grepl("^(median_synthetic|reported_overall_restored):", opes_basis),
        sub("^[^:]+:", "", opes_basis), NA_character_),
      opes_val = coalesce(as.numeric(opes_val_str), stat),
      opes_val_clamped = pmin(pmax(opes_val, -0.995), 0.995),
      yi_opes = atanh(opes_val_clamped),
      vi_opes = case_when(
        grepl("^(median_synthetic|reported_overall_restored):", opes_basis) ~
          ifelse(!is.na(nc) & nc > 3, 1 / (nc - 3), vi),
        TRUE ~ vi
      )
    )
  mod <- rma(yi_opes, vi_opes, data = data, method = "REML")
  pooled_kappa <- tanh(mod$b[1])
  ci_kappa <- tanh(c(mod$ci.lb, mod$ci.ub))
  pi <- predict(mod)
  pi_kappa <- tanh(c(pi$pi.lb, pi$pi.ub))

  cat(sprintf("\n--- %s (OPES, k=%d) ---\n", label, nrow(data)))
  cat(sprintf("Pooled: %.3f (%.3f, %.3f), PI: %.3f to %.3f\n",
              pooled_kappa, ci_kappa[1], ci_kappa[2], pi_kappa[1], pi_kappa[2]))
  cat(sprintf("tau2=%.4f, I2=%.1f%%, Q(df=%d)=%.1f, p=%s\n",
              mod$tau2, mod$I2, mod$k - 1, mod$QE, format.pval(mod$QEp, digits = 3)))
  list(mod = mod, pooled_kappa = pooled_kappa, ci_kappa = ci_kappa, pi_kappa = pi_kappa)
}

run_moderator <- function(data, formula, label) {
  tryCatch({
    mod <- rma.mv(formula, V = data$vi,
                  random = list(~ 1 | dep_id, ~ 1 | comparison_id),
                  data = data, method = "REML")
    cr2 <- coef_test(mod, vcov = "CR2", cluster = data$dep_id)
    cat(sprintf("\n  %s (k=%d):\n", label, nrow(data)))
    for (i in seq_len(nrow(cr2))) {
      cat(sprintf("    %s: beta=%.4f, SE=%.4f, df=%.1f, p=%s\n",
                  rownames(cr2)[i], cr2$beta[i], cr2$SE[i],
                  cr2$df_Satt[i], format.pval(cr2$p_Satt[i], digits = 3)))
    }
    invisible(list(mod = mod, cr2 = cr2))
  }, error = function(e) {
    cat(sprintf("\n  %s: FAILED (%s)\n", label, e$message))
    invisible(NULL)
  })
}

# Summary table helper
summary_rows <- list()
add_row <- function(label, res, k, clust) {
  if (is.null(res)) return()
  summary_rows[[length(summary_rows) + 1]] <<- list(
    label = label, k = k, clust = clust,
    kappa = res$pooled_kappa, ci = res$ci_kappa, pi = res$pi_kappa)
}


# ============================================================
# 3. PRIMARY: KAPPA-ONLY, RVE+CR2
# ============================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PRIMARY ANALYSIS: Kappa-only, RVE+CR2\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

dat_kappa <- dat %>% filter(stat_type == "kappa")
cat(sprintf("\nKappa-only: %d effects\n", nrow(dat_kappa)))

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

validate_opes_selection(dat_kappa)

# --- Inter-reader ---
ir_k <- dat_kappa %>% filter(analytic_stratum == "inter-reader")
ir_res <- run_rve_cr2(ir_k, "Inter-reader PRIMARY (kappa-only)")
add_row("IR kappa RVE", ir_res, ir_res$k, ir_res$n_clusters)

# --- Inter-modality ---
im_k <- dat_kappa %>% filter(analytic_stratum == "inter-modality")
im_res <- run_rve_cr2(im_k, "Inter-modality PRIMARY (kappa-only)")
add_row("IM kappa RVE", im_res, im_res$k, im_res$n_clusters)

# --- OPES ---
opes_k <- dat_kappa %>% filter(opes_include == 1)

ir_opes_k <- opes_k %>% filter(analytic_stratum == "inter-reader")
ir_opes_res <- run_opes_re(ir_opes_k, "Inter-reader kappa")
add_row("IR kappa OPES", ir_opes_res, nrow(ir_opes_k), nrow(ir_opes_k))

im_opes_k <- opes_k %>% filter(analytic_stratum == "inter-modality")
im_opes_res <- run_opes_re(im_opes_k, "Inter-modality kappa")
add_row("IM kappa OPES", im_opes_res, nrow(im_opes_k), nrow(im_opes_k))

intra_opes_k <- opes_k %>% filter(analytic_stratum == "intra-reader")
intra_opes_res <- run_opes_re(intra_opes_k, "Intra-reader kappa")
add_row("Intra kappa OPES", intra_opes_res, nrow(intra_opes_k), nrow(intra_opes_k))

iv_opes_k <- opes_k %>% filter(analytic_stratum == "inter-version")
iv_opes_res <- run_opes_re(iv_opes_k, "Inter-version kappa")
add_row("IV kappa OPES", iv_opes_res, nrow(iv_opes_k), nrow(iv_opes_k))


# ============================================================
# 4. META-REGRESSION: Inter-reader kappa-only
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("META-REGRESSION: Inter-reader kappa-only\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# ---- 4A. UNIVARIABLE (main 5) ----
cat("\n--- Univariable: Main 5 moderators ---\n")

# 1. n_categories (task difficulty — most important)
ir_k_nc <- ir_k %>%
  filter(ncat_group %in% c("2-3", "4-5")) %>%
  mutate(is_full = ifelse(ncat_group == "4-5", 1, 0))
run_moderator(ir_k_nc, yi ~ is_full, "IR.1 n_categories (4-5 vs 2-3)")

# 2. weight_scheme (weighted vs unweighted, classified only)
ir_k_wt <- ir_k %>%
  filter(wt_group %in% c("weighted", "unweighted")) %>%
  mutate(is_weighted = ifelse(wt_group == "weighted", 1, 0))
run_moderator(ir_k_wt, yi ~ is_weighted, "IR.2 weight_scheme (weighted vs UW)")

# 3. Modality (CT vs MRI vs CEUS_US)
ir_k_mod <- ir_k %>%
  filter(mod_group %in% c("CT", "MRI", "CEUS_US")) %>%
  mutate(modality = factor(mod_group, levels = c("CT", "MRI", "CEUS_US")))
run_moderator(ir_k_mod, yi ~ modality, "IR.3 modality (CT/MRI/CEUS_US)")

# 4a. Bosniak version — CT-only subset (PRIMARY: removes modality confounding)
ir_k_ver_ct <- ir_k %>%
  filter(std_modality_1 == "CT",
         std_version_1 %in% c("original", "v2019")) %>%
  mutate(is_v2019 = ifelse(std_version_1 == "v2019", 1, 0))
cat(sprintf("    CT-only version subset: k=%d (original=%d, v2019=%d)\n",
            nrow(ir_k_ver_ct),
            sum(ir_k_ver_ct$is_v2019 == 0),
            sum(ir_k_ver_ct$is_v2019 == 1)))
run_moderator(ir_k_ver_ct, yi ~ is_v2019,
              "IR.4a version CT-only (v2019 vs original)")

# 4b. Bosniak version — all modalities (EXPLORATORY: modality-confounded)
ir_k_ver <- ir_k %>%
  filter(!is.na(ver_group)) %>%
  mutate(is_v2019 = ifelse(ver_group == "v2019", 1, 0))
run_moderator(ir_k_ver, yi ~ is_v2019,
              "IR.4b version all-mod (v2019 vs pre2019, EXPLORATORY)")

# 5. Publication type (fulltext vs abstract)
ir_k_pub <- ir_k %>%
  filter(!is.na(pub_type)) %>%
  mutate(is_abstract = ifelse(pub_type == "abstract", 1, 0))
run_moderator(ir_k_pub, yi ~ is_abstract, "IR.5 pub_type (abstract vs fulltext)")

# ---- 4B. UNIVARIABLE (exploratory) ----
cat("\n--- Univariable: Exploratory ---\n")

# 6. log(n_cysts) — sample size effect
ir_k_lnc <- ir_k %>% filter(!is.na(log_nc))
run_moderator(ir_k_lnc, yi ~ log_nc, "IR.6 log(n_cysts)")

# 7. n_readers (2 vs 3+)
ir_k_nr <- ir_k %>% filter(!is.na(nr_group)) %>%
  mutate(is_3plus = ifelse(nr_group == "3plus", 1, 0))
run_moderator(ir_k_nr, yi ~ is_3plus, "IR.7 n_readers (3+ vs 2)")

# 8. Reader experience (continuous, mean years)
ir_k_exp <- ir_k %>% filter(!is.na(exp_mean))
cat(sprintf("    Experience data available: k=%d/%d\n", nrow(ir_k_exp), nrow(ir_k)))
run_moderator(ir_k_exp, yi ~ exp_mean, "IR.8 experience (continuous)")

# 9. Specialty (subspecialty vs non-subspecialty)
ir_k_spec <- ir_k %>%
  filter(!is.na(spec_group)) %>%
  mutate(is_sub = ifelse(spec_group == "subspecialty", 1, 0))
cat(sprintf("    Specialty data available: k=%d/%d (sub=%d, non=%d)\n",
            nrow(ir_k_spec), nrow(ir_k),
            sum(ir_k_spec$is_sub == 1), sum(ir_k_spec$is_sub == 0)))
run_moderator(ir_k_spec, yi ~ is_sub, "IR.9 specialty (subspecialty vs other)")

# 10. Blinding (blinded vs unclear/not)
ir_k_blind <- ir_k %>% filter(!is.na(is_blinded))
cat(sprintf("    Blinding data available: k=%d/%d (blinded=%d, not=%d)\n",
            nrow(ir_k_blind), nrow(ir_k),
            sum(ir_k_blind$is_blinded == 1), sum(ir_k_blind$is_blinded == 0)))
run_moderator(ir_k_blind, yi ~ is_blinded, "IR.10 blinding (blinded vs other)")

# 11. Reader structure (pooled vs non-pooled)
ir_k_pool <- ir_k %>% filter(!is.na(is_pooled))
cat(sprintf("    Reader structure: k=%d (pooled=%d, non-pooled=%d)\n",
            nrow(ir_k_pool), sum(ir_k_pool$is_pooled == 1),
            sum(ir_k_pool$is_pooled == 0)))
run_moderator(ir_k_pool, yi ~ is_pooled, "IR.11 reader_structure (pooled vs other)")

# 12. Publication year (continuous)
ir_k_yr <- ir_k %>% filter(!is.na(pub_year))
run_moderator(ir_k_yr, yi ~ pub_year, "IR.12 publication year")

# 13. QAREL quality (Low vs Moderate/High)
ir_k_rob <- ir_k %>% filter(!is.na(is_low_rob))
cat(sprintf("    QAREL quality: k=%d (Low=%d, Moderate/High=%d)\n",
            nrow(ir_k_rob), sum(ir_k_rob$is_low_rob == 1),
            sum(ir_k_rob$is_low_rob == 0)))
run_moderator(ir_k_rob, yi ~ is_low_rob, "IR.13 QAREL quality (Low vs Mod/High)")

# 13b. QAREL yes-count (continuous)
ir_k_qs <- ir_k %>% filter(!is.na(qarel_score))
run_moderator(ir_k_qs, yi ~ qarel_score, "IR.13b QAREL yes-count (continuous)")

# 14. Study design (prospective vs retrospective)
ir_k_des <- ir_k %>% filter(!is.na(is_prospective))
cat(sprintf("    IR design data available: k=%d (retro=%d, prosp=%d)\n",
            nrow(ir_k_des), sum(ir_k_des$is_prospective == 0),
            sum(ir_k_des$is_prospective == 1)))
run_moderator(ir_k_des, yi ~ is_prospective, "IR.14 design (prospective vs retrospective)")

# 15. Region (North America / Asia / Europe)
ir_k_reg <- ir_k %>%
  filter(!is.na(region_cat)) %>%
  mutate(region_cat = factor(region_cat, levels = c("North_America", "Asia", "Europe")))
cat(sprintf("    IR region data available: k=%d (NA=%d, Asia=%d, Europe=%d)\n",
            nrow(ir_k_reg),
            sum(ir_k_reg$region_cat == "North_America"),
            sum(ir_k_reg$region_cat == "Asia"),
            sum(ir_k_reg$region_cat == "Europe")))
run_moderator(ir_k_reg, yi ~ region_cat, "IR.15 region (NA/Asia/Europe)")

# ---- 4C. SENSITIVITY: weight 3-level (weighted/UW/unspecified) ----
cat("\n--- Sensitivity: weight_scheme 3-level ---\n")
ir_k_wt3 <- ir_k %>%
  filter(!is.na(wt_3level)) %>%
  mutate(wt3 = factor(wt_3level, levels = c("unweighted", "weighted", "unspecified")))
run_moderator(ir_k_wt3, yi ~ wt3, "IR.S weight_scheme 3-level")

# ---- 4D. MULTIVARIABLE (3 models) ----
cat("\n--- Multivariable meta-regression ---\n")

# Model 1a: n_categories + version — CT-only (PRIMARY, no modality confound)
ir_k_mv1a <- ir_k %>%
  filter(std_modality_1 == "CT",
         std_version_1 %in% c("original", "v2019"),
         ncat_group %in% c("2-3", "4-5")) %>%
  mutate(is_full = ifelse(ncat_group == "4-5", 1, 0),
         is_v2019 = ifelse(std_version_1 == "v2019", 1, 0))
run_moderator(ir_k_mv1a, yi ~ is_full + is_v2019,
              "IR.MV1a CT-only: ncat + version")

# Model 1b: n_categories + modality + version — all data (EXPLORATORY)
ir_k_mv1b <- ir_k %>%
  filter(ncat_group %in% c("2-3", "4-5"),
         mod_group %in% c("CT", "MRI", "CEUS_US"),
         !is.na(ver_group)) %>%
  mutate(is_full = ifelse(ncat_group == "4-5", 1, 0),
         modality = factor(mod_group, levels = c("CT", "MRI", "CEUS_US")),
         is_v2019 = ifelse(ver_group == "v2019", 1, 0))
run_moderator(ir_k_mv1b, yi ~ is_full + modality + is_v2019,
              "IR.MV1b all-mod: ncat + modality + version (EXPLORATORY)")

# Model 2: n_categories + weight_scheme (methodological)
ir_k_mv2 <- ir_k %>%
  filter(ncat_group %in% c("2-3", "4-5"),
         wt_group %in% c("weighted", "unweighted")) %>%
  mutate(is_full = ifelse(ncat_group == "4-5", 1, 0),
         is_weighted = ifelse(wt_group == "weighted", 1, 0))
run_moderator(ir_k_mv2, yi ~ is_full + is_weighted,
              "IR.MV2 ncat + weight_scheme")

# Model 3: experience + n_categories (clinical + task difficulty)
ir_k_mv3 <- ir_k %>%
  filter(!is.na(exp_mean), ncat_group %in% c("2-3", "4-5")) %>%
  mutate(is_full = ifelse(ncat_group == "4-5", 1, 0))
cat(sprintf("    MV3 experience + ncat: k=%d\n", nrow(ir_k_mv3)))
run_moderator(ir_k_mv3, yi ~ exp_mean + is_full,
              "IR.MV3 experience + ncat")


# ============================================================
# 4E. META-REGRESSION: Inter-modality kappa-only
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("META-REGRESSION: Inter-modality kappa-only\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\n--- Univariable ---\n")

# 1. n_categories
im_k_nc <- im_k %>%
  filter(ncat_group %in% c("2-3", "4-5")) %>%
  mutate(is_full = ifelse(ncat_group == "4-5", 1, 0))
run_moderator(im_k_nc, yi ~ is_full, "IM.1 n_categories (4-5 vs 2-3)")

# NOTE: version moderator REMOVED for inter-modality.
# Original Bosniak = CT-only; v2019 = CT+MRI. In IM comparisons,
# "version" is inseparable from classification scope (which modalities
# the schema was designed to cover). CT-CEUS comparisons further
# complicate interpretation since neither version covers CEUS.

# 2. modality pair (CT_MRI vs CT_CEUS vs other)
im_k_mp <- im_k %>%
  filter(mod_pair %in% c("CT_MRI", "CT_CEUS", "other")) %>%
  mutate(pair = factor(mod_pair, levels = c("CT_MRI", "CT_CEUS", "other")))
run_moderator(im_k_mp, yi ~ pair, "IM.2 modality_pair (CT-MRI/CT-CEUS/other)")

# 3. log(n_cysts)
im_k_lnc <- im_k %>% filter(!is.na(log_nc))
run_moderator(im_k_lnc, yi ~ log_nc, "IM.3 log(n_cysts)")

# 4. publication type
im_k_pub <- im_k %>%
  filter(!is.na(pub_type)) %>%
  mutate(is_abstract = ifelse(pub_type == "abstract", 1, 0))
run_moderator(im_k_pub, yi ~ is_abstract, "IM.4 pub_type (abstract vs fulltext)")
cat("    Interpretation: not interpretable because Satterthwaite df is below 4 (abstract level: two studies and two dependency clusters).\n")

# 5. Reader experience (continuous)
im_k_exp <- im_k %>% filter(!is.na(exp_mean))
cat(sprintf("    IM experience data available: k=%d/%d\n", nrow(im_k_exp), nrow(im_k)))
run_moderator(im_k_exp, yi ~ exp_mean, "IM.5 experience (continuous)")

# 6. Specialty (subspecialty vs other)
im_k_spec <- im_k %>%
  filter(!is.na(spec_group)) %>%
  mutate(is_sub = ifelse(spec_group == "subspecialty", 1, 0))
run_moderator(im_k_spec, yi ~ is_sub, "IM.6 specialty (subspecialty vs other)")

# 7. Blinding (blinded vs other)
im_k_blind <- im_k %>% filter(!is.na(is_blinded))
run_moderator(im_k_blind, yi ~ is_blinded, "IM.7 blinding (blinded vs other)")

# 8. Publication year (continuous)
im_k_yr <- im_k %>% filter(!is.na(pub_year))
run_moderator(im_k_yr, yi ~ pub_year, "IM.8 publication year")

# 9. Weight scheme (weighted vs unweighted)
im_k_wt <- im_k %>%
  filter(wt_group %in% c("weighted", "unweighted")) %>%
  mutate(is_weighted = ifelse(wt_group == "weighted", 1, 0))
run_moderator(im_k_wt, yi ~ is_weighted, "IM.9 weight_scheme (weighted vs UW)")

# 10. n_readers (2 vs 3+ — for IM this may vary with 1-reader designs)
im_k_nr <- im_k %>%
  filter(!is.na(n_readers)) %>%
  mutate(nr_cat = case_when(
    n_readers == 1 ~ "1",
    n_readers == 2 ~ "2",
    n_readers >= 3 ~ "3plus",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(nr_cat))
if (length(unique(im_k_nr$nr_cat)) >= 2) {
  im_k_nr$nr_cat <- factor(im_k_nr$nr_cat, levels = c("1", "2", "3plus"))
  run_moderator(im_k_nr, yi ~ nr_cat, "IM.10 n_readers (1/2/3+)")
}

# 11. QAREL quality (Low vs Moderate/High)
im_k_rob <- im_k %>% filter(!is.na(is_low_rob))
cat(sprintf("    IM QAREL quality: k=%d (Low=%d, Moderate/High=%d)\n",
            nrow(im_k_rob), sum(im_k_rob$is_low_rob == 1),
            sum(im_k_rob$is_low_rob == 0)))
run_moderator(im_k_rob, yi ~ is_low_rob, "IM.11 QAREL quality (Low vs Mod/High)")

# 11b. QAREL yes-count (continuous)
im_k_qs <- im_k %>% filter(!is.na(qarel_score))
run_moderator(im_k_qs, yi ~ qarel_score, "IM.11b QAREL yes-count (continuous)")

# 12. Study design (prospective vs retrospective)
im_k_des <- im_k %>% filter(!is.na(is_prospective))
cat(sprintf("    IM design data available: k=%d (retro=%d, prosp=%d)\n",
            nrow(im_k_des), sum(im_k_des$is_prospective == 0),
            sum(im_k_des$is_prospective == 1)))
run_moderator(im_k_des, yi ~ is_prospective, "IM.12 design (prospective vs retrospective)")

# 13. Region (North America / Asia / Europe)
im_k_reg <- im_k %>%
  filter(!is.na(region_cat)) %>%
  mutate(region_cat = factor(region_cat, levels = c("North_America", "Asia", "Europe")))
cat(sprintf("    IM region data available: k=%d (NA=%d, Asia=%d, Europe=%d)\n",
            nrow(im_k_reg),
            sum(im_k_reg$region_cat == "North_America"),
            sum(im_k_reg$region_cat == "Asia"),
            sum(im_k_reg$region_cat == "Europe")))
run_moderator(im_k_reg, yi ~ region_cat, "IM.13 region (NA/Asia/Europe)")


# ============================================================
# 5. SENSITIVITY 1: Weighted kappa only (+ AT-derived LW)
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SENSITIVITY 1: Weighted kappa only (+ AT-derived LW)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Prepare AT-derived LW dataset through the same pipeline as primary
dat_at_lw <- comp %>%
  filter(exclude == "at_derived_lw_sensitivity") %>%
  left_join(stud %>% select(study_id, year, publication_type, study_design),
            by = "study_id") %>%
  left_join(qarel %>% select(study_id, overall_rob, qarel_yes_count),
            by = "study_id") %>%
  mutate(
    kappa = as.numeric(kappa),
    nc = as.numeric(n_cysts),
    stat_type = "kappa",
    stat = kappa,
    stat_clamped = pmin(pmax(stat, -0.995), 0.995),
    yi = atanh(stat_clamped),
    vi_n = ifelse(!is.na(nc) & nc > 3, 1 / (nc - 3), NA_real_),
    vi = vi_n,
    wt_group = "weighted"
  ) %>%
  filter(!is.na(yi) & !is.na(vi) & vi > 0)
cat(sprintf("AT-derived LW rows loaded: %d\n", nrow(dat_at_lw)))

at_lw_ir <- dat_at_lw %>% filter(analytic_stratum == "inter-reader")
at_lw_im <- dat_at_lw %>% filter(analytic_stratum == "inter-modality")

# Weighted sensitivity: author-reported weighted + AT-derived LW
ir_wk <- bind_rows(
  ir_k %>% filter(wt_group == "weighted"),
  at_lw_ir
)
ir_wk_res <- run_rve_cr2(ir_wk, "IR weighted kappa (+AT-derived LW)")
add_row("IR weighted RVE", ir_wk_res, ir_wk_res$k, ir_wk_res$n_clusters)

im_wk <- bind_rows(
  im_k %>% filter(wt_group == "weighted"),
  at_lw_im
)
im_wk_res <- run_rve_cr2(im_wk, "IM weighted kappa (+AT-derived LW)")


# ============================================================
# 6. SENSITIVITY 2: Unweighted kappa only
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SENSITIVITY 2: Unweighted kappa only\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

ir_uk <- ir_k %>% filter(wt_group == "unweighted")
ir_uk_res <- run_rve_cr2(ir_uk, "IR unweighted kappa")
add_row("IR unweighted RVE", ir_uk_res, ir_uk_res$k, ir_uk_res$n_clusters)

im_uk <- im_k %>% filter(wt_group == "unweighted")
im_uk_res <- run_rve_cr2(im_uk, "IM unweighted kappa")
add_row("IM unweighted RVE", im_uk_res, im_uk_res$k, im_uk_res$n_clusters)


# ============================================================
# 7. SENSITIVITY 3: 5-category kappa only
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SENSITIVITY 3: 5-category kappa only\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

ir_5c <- ir_k %>% filter(n_categories == 5)
ir_5c_res <- run_rve_cr2(ir_5c, "IR 5-cat kappa")
add_row("IR 5-cat RVE", ir_5c_res, ir_5c_res$k, ir_5c_res$n_clusters)

im_5c <- im_k %>% filter(n_categories == 5)
im_5c_res <- run_rve_cr2(im_5c, "IM 5-cat kappa")
add_row("IM 5-cat RVE", im_5c_res, im_5c_res$k, im_5c_res$n_clusters)


# ============================================================
# 8. EXPLORATORY: All metrics pooled
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("EXPLORATORY: All metrics pooled (kappa + Gwet + ICC)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

ir_all <- dat %>% filter(analytic_stratum == "inter-reader")
ir_all_res <- run_rve_cr2(ir_all, "IR all metrics")
add_row("IR all-metric RVE", ir_all_res, ir_all_res$k, ir_all_res$n_clusters)

# stat_type moderator
ir_all_st <- ir_all %>%
  mutate(stat_type = factor(stat_type, levels = c("kappa", "gwet", "icc")))
run_moderator(ir_all_st, yi ~ stat_type, "stat_type (kappa/gwet/icc)")

im_all <- dat %>% filter(analytic_stratum == "inter-modality")
im_all_res <- run_rve_cr2(im_all, "IM all metrics")
add_row("IM all-metric RVE", im_all_res, im_all_res$k, im_all_res$n_clusters)


# ============================================================
# 9. EXPLORATORY: Small strata all-effects (kappa-only)
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("EXPLORATORY: Small strata all-effects (kappa-only)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

intra_all_k <- dat_kappa %>% filter(analytic_stratum == "intra-reader")
intra_rve <- run_rve_cr2(intra_all_k, "Intra-reader kappa all-effects")

iv_all_k <- dat_kappa %>% filter(analytic_stratum == "inter-version")
iv_rve <- run_rve_cr2(iv_all_k, "Inter-version kappa all-effects")


# ============================================================
# 10. RVE vs OPES COMPARISON
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("RVE vs OPES COMPARISON (kappa-only)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

compare <- function(rve, opes, label) {
  if (is.null(rve) || is.null(opes)) return()
  diff <- (opes$pooled_kappa - rve$pooled_kappa) / rve$pooled_kappa * 100
  overlap <- !(opes$ci_kappa[2] < rve$ci_kappa[1] | rve$ci_kappa[2] < opes$ci_kappa[1])
  cat(sprintf("\n%s:\n  RVE: %.3f (%.3f, %.3f)\n  OPES: %.3f (%.3f, %.3f)\n  Diff: %+.1f%%, CI overlap: %s\n",
              label,
              rve$pooled_kappa, rve$ci_kappa[1], rve$ci_kappa[2],
              opes$pooled_kappa, opes$ci_kappa[1], opes$ci_kappa[2],
              diff, ifelse(overlap, "YES", "NO")))
}

compare(ir_res, ir_opes_res, "Inter-reader")
compare(im_res, im_opes_res, "Inter-modality")


# ============================================================
# 11. SUMMARY TABLE
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SUMMARY TABLE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat(sprintf("\n%-25s %5s %6s %8s %18s %18s\n",
            "Analysis", "k", "clust", "est", "95% CI", "95% PI"))
cat(paste(rep("-", 85), collapse = ""), "\n")

for (r in summary_rows) {
  cat(sprintf("%-25s %5d %6d %8.3f  [%.3f, %.3f]  [%.3f, %.3f]\n",
              r$label, r$k, r$clust, r$kappa,
              r$ci[1], r$ci[2], r$pi[1], r$pi[2]))
}

cat("\nDone.\n")
