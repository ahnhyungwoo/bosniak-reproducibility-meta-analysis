## ============================================================
## Bosniak Classification Reproducibility SRMA
## Main Meta-Analysis â€” Revised Framework
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

# Source verification during revision showed that Chang 2015 was a conference
# abstract rather than a full journal article. Apply the corrected publication
# type without changing the deposited source table.
stud <- stud %>%
  mutate(
    publication_type = if_else(
      study_id == "Chang_2015",
      "conference_abstract",
      publication_type
    )
  )

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

    # Region (3-level: North_America, Asia, Europe; Other/empty â†’ NA)
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

# 1. n_categories (task difficulty â€” most important)×_6¶‰žËkºwµç@€€€€€€¹É½Ü¡¥É}­}Á½½°¤°ÍÕ´¡¥É}­}Á½½°‘¥Í}Á½½±•€ôô€Ä¤°4(€€€€€€€€€€€ÍÕ´¡¥É}­}Á½½°‘¥Í}Á½½±•€ôô€À¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}Á½½°°å¤ø¥Í}Á½½±•°€‰%H¸ÄÄÉ•…‘•É}ÍÑÉÕÑÕÉ”€¡Á½½±•ÙÌ½Ñ¡•È¤ˆ¤4(4(Œ€ÄÈ¸AÕ‰±¥…Ñ¥½¸å•…È€¡½¹Ñ¥¹Õ½ÕÌ¤4)¥É}­}åÈ€ð´¥É}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡ÁÕ‰}å•…È¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}åÈ°å¤øÁÕ‰}å•…È°€‰%H¸ÄÈÁÕ‰±¥…Ñ¥½¸å•…Èˆ¤4(4(Œ€ÄÌ¸EI0ÅÕ…±¥Ñä€¡1½ÜÙÌ5½‘•É…Ñ”½!¥ ¤4)¥É}­}É½ˆ€ð´¥É}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡¥Í}±½Ý}É½ˆ¤¤4)…Ð¡ÍÁÉ¥¹Ñ˜ ˆ€€€EI0ÅÕ…±¥Ñäè¬ô•€¡1½Üô•°5½‘•É…Ñ”½!¥ ô•¥q¸ˆ°4(€€€€€€€€€€€¹É½Ü¡¥É}­}É½ˆ¤°ÍÕ´¡¥É}­}É½ˆ‘¥Í}±½Ý}É½ˆ€ôô€Ä¤°4(€€€€€€€€€€€ÍÕ´¡¥É}­}É½ˆ‘¥Í}±½Ý}É½ˆ€ôô€À¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}É½ˆ°å¤ø¥Í}±½Ý}É½ˆ°€‰%H¸ÄÌEI0ÅÕ…±¥Ñä€¡1½ÜÙÌ5½½!¥ ¤ˆ¤4(4(Œ€ÄÍˆ¸EI0å•Ìµ½Õ¹Ð€¡½¹Ñ¥¹Õ½ÕÌ¤4)¥É}­}ÅÌ€ð´¥É}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡Å…É•±}Í½É”¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}ÅÌ°å¤øÅ…É•±}Í½É”°€‰%H¸ÄÍˆEI0å•Ìµ½Õ¹Ð€¡½¹Ñ¥¹Õ½ÕÌ¤ˆ¤4(4(Œ€ÄÐ¸MÑÕ‘ä‘•Í¥¸€¡ÁÉ½ÍÁ•Ñ¥Ù”ÙÌÉ•ÑÉ½ÍÁ•Ñ¥Ù”¤4)¥É}­}‘•Ì€ð´¥É}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡¥Í}ÁÉ½ÍÁ•Ñ¥Ù”¤¤4)…Ð¡ÍÁÉ¥¹Ñ˜ ˆ€€€%H‘•Í¥¸‘…Ñ„…Ù…¥±…‰±”è¬ô•€¡É•ÑÉ¼ô•°ÁÉ½ÍÀô•¥q¸ˆ°4(€€€€€€€€€€€¹É½Ü¡¥É}­}‘•Ì¤°ÍÕ´¡¥É}­}‘•Ì‘¥Í}ÁÉ½ÍÁ•Ñ¥Ù”€ôô€À¤°4(€€€€€€€€€€€ÍÕ´¡¥É}­}‘•Ì‘¥Í}ÁÉ½ÍÁ•Ñ¥Ù”€ôô€Ä¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}‘•Ì°å¤ø¥Í}ÁÉ½ÍÁ•Ñ¥Ù”°€‰%H¸ÄÐ‘•Í¥¸€¡ÁÉ½ÍÁ•Ñ¥Ù”ÙÌÉ•ÑÉ½ÍÁ•Ñ¥Ù”¤ˆ¤4(4(Œ€ÄÔ¸I•¥½¸€¡9½ÉÑ µ•É¥„€¼Í¥„€¼ÕÉ½Á”¤4)¥É}­}É•œ€ð´¥É}¬€”ø”4(€™¥±Ñ•È …¥Ì¹¹„¡É•¥½¹}…Ð¤¤€”ø”4(€µÕÑ…Ñ”¡É•¥½¹}…Ð€ô™…Ñ½È¡É•¥½¹}…Ð°±•Ù•±Ì€ôŒ ‰9½ÉÑ¡}µ•É¥„ˆ°€‰Í¥„ˆ°€‰ÕÉ½Á”ˆ¤¤¤4)…Ð¡ÍÁÉ¥¹Ñ˜ ˆ€€€%HÉ•¥½¸‘…Ñ„…Ù…¥±…‰±”è¬ô•€¡9ô•°Í¥„ô•°ÕÉ½Á”ô•¥q¸ˆ°4(€€€€€€€€€€€¹É½Ü¡¥É}­}É•œ¤°4(€€€€€€€€€€€ÍÕ´¡¥É}­}É•œ‘É•¥½¹}…Ð€ôô€‰9½ÉÑ¡}µ•É¥„ˆ¤°4(€€€€€€€€€€€ÍÕ´¡¥É}­}É•œ‘É•¥½¹}…Ð€ôô€‰Í¥„ˆ¤°4(€€€€€€€€€€€ÍÕ´¡¥É}­}É•œ‘É•¥½¹}…Ð€ôô€‰ÕÉ½Á”ˆ¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}É•œ°å¤øÉ•¥½¹}…Ð°€‰%H¸ÄÔÉ•¥½¸€¡9½Í¥„½ÕÉ½Á”¤ˆ¤4(4(Œ€´´´´€Ñ¸M9M%Q%Y%QdèÝ•¥¡Ð€Ìµ±•Ù•°€¡Ý•¥¡Ñ•½U\½Õ¹ÍÁ•¥™¥•¤€´´´´4)…Ð ‰q¸´´´M•¹Í¥Ñ¥Ù¥ÑäèÝ•¥¡Ñ}Í¡•µ”€Ìµ±•Ù•°€´´µq¸ˆ¤4)¥É}­}ÝÐÌ€ð´¥É}¬€”ø”4(€™¥±Ñ•È …¥Ì¹¹„¡ÝÑ|Í±•Ù•°¤¤€”ø”4(€µÕÑ…Ñ”¡ÝÐÌ€ô™…Ñ½È¡ÝÑ|Í±•Ù•°°±•Ù•±Ì€ôŒ ‰Õ¹Ý•¥¡Ñ•ˆ°€‰Ý•¥¡Ñ•ˆ°€‰Õ¹ÍÁ•¥™¥•ˆ¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}ÝÐÌ°å¤øÝÐÌ°€‰%H¹LÝ•¥¡Ñ}Í¡•µ”€Ìµ±•Ù•°ˆ¤4(4(Œ€´´´´€Ñ¸5U1Q%YI%	1€ Ìµ½‘•±Ì¤€´´´´4)…Ð ‰q¸´´´5Õ±Ñ¥Ù…É¥…‰±”µ•Ñ„µÉ•É•ÍÍ¥½¸€´´µq¸ˆ¤4(4(Œ5½‘•°€Å„è¹}…Ñ•½É¥•Ì€¬Ù•ÉÍ¥½¸ƒŠPPµ½¹±ä€¡AI%5Id°¹¼µ½‘…±¥Ñä½¹™½Õ¹¤4)¥É}­}µØÅ„€ð´¥É}¬€”ø”4(€™¥±Ñ•È¡ÍÑ‘}µ½‘…±¥Ñå|Ä€ôô€‰Pˆ°4(€€€€€€€€ÍÑ‘}Ù•ÉÍ¥½¹|Ä€•¥¸”Œ ‰½É¥¥¹…°ˆ°€‰ØÈÀÄäˆ¤°4(€€€€€€€€¹…Ñ}É½ÕÀ€•¥¸”Œ ˆÈ´Ìˆ°€ˆÐ´Ôˆ¤¤€”ø”4(€µÕÑ…Ñ”¡¥Í}™Õ±°€ô¥™•±Í”¡¹…Ñ}É½ÕÀ€ôô€ˆÐ´Ôˆ°€Ä°€À¤°4(€€€€€€€€¥Í}ØÈÀÄä€ô¥™•±Í”¡ÍÑ‘}Ù•ÉÍ¥½¹|Ä€ôô€‰ØÈÀÄäˆ°€Ä°€À¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}µØÅ„°å¤ø¥Í}™Õ±°€¬¥Í}ØÈÀÄä°4(€€€€€€€€€€€€€€‰%H¹5XÅ„Pµ½¹±äè¹…Ð€¬Ù•ÉÍ¥½¸ˆ¤4(4(Œ5½‘•°€Åˆè¹}…Ñ•½É¥•Ì€¬µ½‘…±¥Ñä€¬Ù•ÉÍ¥½¸ƒŠP…±°‘…Ñ„€¡aA1=IQ=Id¤4)¥É}­}µØÅˆ€ð´¥É}¬€”ø”4(€™¥±Ñ•È¡¹…Ñ}É½ÕÀ€•¥¸”Œ ˆÈ´Ìˆ°€ˆÐ´Ôˆ¤°4(€€€€€€€€µ½‘}É½ÕÀ€•¥¸”Œ ‰Pˆ°€‰5I$ˆ°€‰UM}ULˆ¤°4(€€€€€€€€€…¥Ì¹¹„¡Ù•É}É½ÕÀ¤¤€”ø”4(€µÕÑ…Ñ”¡¥Í}™Õ±°€ô¥™•±Í”¡¹…Ñ}É½ÕÀ€ôô€ˆÐ´Ôˆ°€Ä°€À¤°4(€€€€€€€€µ½‘…±¥Ñä€ô™…Ñ½È¡µ½‘}É½ÕÀ°±•Ù•±Ì€ôŒ ‰Pˆ°€‰5I$ˆ°€‰UM}ULˆ¤¤°4(€€€€€€€€¥Í}ØÈÀÄä€ô¥™•±Í”¡Ù•É}É½ÕÀ€ôô€‰ØÈÀÄäˆ°€Ä°€À¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}µØÅˆ°å¤ø¥Í}™Õ±°€¬µ½‘…±¥Ñä€¬¥Í}ØÈÀÄä°4(€€€€€€€€€€€€€€‰%H¹5XÅˆ…±°µµ½è¹…Ð€¬µ½‘…±¥Ñä€¬Ù•ÉÍ¥½¸€¡aA1=IQ=Id¤ˆ¤4(4(Œ5½‘•°€Èè¹}…Ñ•½É¥•Ì€¬Ý•¥¡Ñ}Í¡•µ”€¡µ•Ñ¡½‘½±½¥…°¤4)¥É}­}µØÈ€ð´¥É}¬€”ø”4(€™¥±Ñ•È¡¹…Ñ}É½ÕÀ€•¥¸”Œ ˆÈ´Ìˆ°€ˆÐ´Ôˆ¤°4(€€€€€€€€ÝÑ}É½ÕÀ€•¥¸”Œ ‰Ý•¥¡Ñ•ˆ°€‰Õ¹Ý•¥¡Ñ•ˆ¤¤€”ø”4(€µÕÑ…Ñ”¡¥Í}™Õ±°€ô¥™•±Í”¡¹…Ñ}É½ÕÀ€ôô€ˆÐ´Ôˆ°€Ä°€À¤°4(€€€€€€€€¥Í}Ý•¥¡Ñ•€ô¥™•±Í”¡ÝÑ}É½ÕÀ€ôô€‰Ý•¥¡Ñ•ˆ°€Ä°€À¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}µØÈ°å¤ø¥Í}™Õ±°€¬¥Í}Ý•¥¡Ñ•°4(€€€€€€€€€€€€€€‰%H¹5XÈ¹…Ð€¬Ý•¥¡Ñ}Í¡•µ”ˆ¤4(4(Œ5½‘•°€Ìè•áÁ•É¥•¹”€¬¹}…Ñ•½É¥•Ì€¡±¥¹¥…°€¬Ñ…Í¬‘¥™™¥Õ±Ñä¤4)¥É}­}µØÌ€ð´¥É}¬€”ø”4(€™¥±Ñ•È …¥Ì¹¹„¡•áÁ}µ•…¸¤°¹…Ñ}É½ÕÀ€•¥¸”Œ ˆÈ´Ìˆ°€ˆÐ´Ôˆ¤¤€”ø”4(€µÕÑ…Ñ”¡¥Í}™Õ±°€ô¥™•±Í”¡¹…Ñ}É½ÕÀ€ôô€ˆÐ´Ôˆ°€Ä°€À¤¤4)…Ð¡ÍÁÉ¥¹Ñ˜ ˆ€€€5XÌ•áÁ•É¥•¹”€¬¹…Ðè¬ô•‘q¸ˆ°¹É½Ü¡¥É}­}µØÌ¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}­}µØÌ°å¤ø•áÁ}µ•…¸€¬¥Í}™Õ±°°4(€€€€€€€€€€€€€€‰%H¹5XÌ•áÁ•É¥•¹”€¬¹…Ðˆ¤4(4(4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(Œ€Ñ¸5QµIIMM%=8è%¹Ñ•Èµµ½‘…±¥Ñä­…ÁÁ„µ½¹±ä4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(4)…Ð ‰q¸ˆ°Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4)…Ð ‰5QµIIMM%=8è%¹Ñ•Èµµ½‘…±¥Ñä­…ÁÁ„µ½¹±åq¸ˆ¤4)…Ð¡Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4(4)…Ð ‰q¸´´´U¹¥Ù…É¥…‰±”€´´µq¸ˆ¤4(4(Œ€Ä¸¹}…Ñ•½É¥•Ì4)¥µ}­}¹Œ€ð´¥µ}¬€”ø”4(€™¥±Ñ•È¡¹…Ñ}É½ÕÀ€•¥¸”Œ ˆÈ´Ìˆ°€ˆÐ´Ôˆ¤¤€”ø”4(€µÕÑ…Ñ”¡¥Í}™Õ±°€ô¥™•±Í”¡¹…Ñ}É½ÕÀ€ôô€ˆÐ´Ôˆ°€Ä°€À¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}¹Œ°å¤ø¥Í}™Õ±°°€‰%4¸Ä¹}…Ñ•½É¥•Ì€ Ð´ÔÙÌ€È´Ì¤ˆ¤4(4(Œ9=QèÙ•ÉÍ¥½¸µ½‘•É…Ñ½ÈI5=Y™½È¥¹Ñ•Èµµ½‘…±¥Ñä¸4(Œ=É¥¥¹…°	½Í¹¥…¬€ôPµ½¹±äìØÈÀÄä€ôP­5I$¸%¸%4½µÁ…É¥Í½¹Ì°4(Œ€‰Ù•ÉÍ¥½¸ˆ¥Ì¥¹Í•Á…É…‰±”™É½´±…ÍÍ¥™¥…Ñ¥½¸Í½Á”€¡Ý¡¥ µ½‘…±¥Ñ¥•Ì4(ŒÑ¡”Í¡•µ„Ý…Ì‘•Í¥¹•Ñ¼½Ù•È¤¸PµUL½µÁ…É¥Í½¹Ì™ÕÉÑ¡•È4(Œ½µÁ±¥…Ñ”¥¹Ñ•ÉÁÉ•Ñ…Ñ¥½¸Í¥¹”¹•¥Ñ¡•ÈÙ•ÉÍ¥½¸½Ù•ÉÌUL¸4(4(Œ€È¸µ½‘…±¥ÑäÁ…¥È€¡Q}5I$ÙÌQ}ULÙÌ½Ñ¡•È¤4)¥µ}­}µÀ€ð´¥µ}¬€”ø”4(€™¥±Ñ•È¡µ½‘}Á…¥È€•¥¸”Œ ‰Q}5I$ˆ°€‰Q}ULˆ°€‰½Ñ¡•Èˆ¤¤€”ø”4(€µÕÑ…Ñ”¡Á…¥È€ô™…Ñ½È¡µ½‘}Á…¥È°±•Ù•±Ì€ôŒ ‰Q}5I$ˆ°€‰Q}ULˆ°€‰½Ñ¡•Èˆ¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}µÀ°å¤øÁ…¥È°€‰%4¸Èµ½‘…±¥Ñå}Á…¥È€¡Pµ5I$½PµUL½½Ñ¡•È¤ˆ¤4(4(Œ€Ì¸±½œ¡¹}åÍÑÌ¤4)¥µ}­}±¹Œ€ð´¥µ}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡±½}¹Œ¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}±¹Œ°å¤ø±½}¹Œ°€‰%4¸Ì±½œ¡¹}åÍÑÌ¤ˆ¤4(4(Œ€Ð¸ÁÕ‰±¥…Ñ¥½¸ÑåÁ”4)¥µ}­}ÁÕˆ€ð´¥µ}¬€”ø”(€™¥±Ñ•È …¥Ì¹¹„¡ÁÕ‰}ÑåÁ”¤¤€”ø”(€µÕÑ…Ñ”¡¥Í}…‰ÍÑÉ…Ð€ô¥™•±Í”¡ÁÕ‰}ÑåÁ”€ôô€‰…‰ÍÑÉ…Ðˆ°€Ä°€À¤¤)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}ÁÕˆ°å¤ø¥Í}…‰ÍÑÉ…Ð°€‰%4¸ÐÁÕ‰}ÑåÁ”€¡…‰ÍÑÉ…ÐÙÌ™Õ±±Ñ•áÐ¤ˆ¤)…Ð ˆ€€€%¹Ñ•ÉÁÉ•Ñ…Ñ¥½¸è¹½Ð¥¹Ñ•ÉÁÉ•Ñ…‰±”‰•…ÕÍ”Ñ¡”…‰ÍÑÉ…Ð±•Ù•°½¹Ñ…¥¹Ì½¹±ä½¹”ÍÑÕ‘ä…¹½¹”‘•Á•¹‘•¹ä±ÕÍÑ•È¹q¸ˆ¤(4(Œ€Ô¸I•…‘•È•áÁ•É¥•¹”€¡½¹Ñ¥¹Õ½ÕÌ¤4)¥µ}­}•áÀ€ð´¥µ}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡•áÁ}µ•…¸¤¤4)…Ð¡ÍÁÉ¥¹Ñ˜ ˆ€€€%4•áÁ•É¥•¹”‘…Ñ„…Ù…¥±…‰±”è¬ô•¼•‘q¸ˆ°¹É½Ü¡¥µ}­}•áÀ¤°¹É½Ü¡¥µ}¬¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}•áÀ°å¤ø•áÁ}µ•…¸°€‰%4¸Ô•áÁ•É¥•¹”€¡½¹Ñ¥¹Õ½ÕÌ¤ˆ¤4(4(Œ€Ø¸MÁ•¥…±Ñä€¡ÍÕ‰ÍÁ•¥…±ÑäÙÌ½Ñ¡•È¤4)¥µ}­}ÍÁ•Œ€ð´¥µ}¬€”ø”4(€™¥±Ñ•È …¥Ì¹¹„¡ÍÁ•}É½ÕÀ¤¤€”ø”4(€µÕÑ…Ñ”¡¥Í}ÍÕˆ€ô¥™•±Í”¡ÍÁ•}É½ÕÀ€ôô€‰ÍÕ‰ÍÁ•¥…±Ñäˆ°€Ä°€À¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}ÍÁ•Œ°å¤ø¥Í}ÍÕˆ°€‰%4¸ØÍÁ•¥…±Ñä€¡ÍÕ‰ÍÁ•¥…±ÑäÙÌ½Ñ¡•È¤ˆ¤4(4(Œ€Ü¸	±¥¹‘¥¹œ€¡‰±¥¹‘•ÙÌ½Ñ¡•È¤4)¥µ}­}‰±¥¹€ð´¥µ}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡¥Í}‰±¥¹‘•¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}‰±¥¹°å¤ø¥Í}‰±¥¹‘•°€‰%4¸Ü‰±¥¹‘¥¹œ€¡‰±¥¹‘•ÙÌ½Ñ¡•È¤ˆ¤4(4(Œ€à¸AÕ‰±¥…Ñ¥½¸å•…È€¡½¹Ñ¥¹Õ½ÕÌ¤4)¥µ}­}åÈ€ð´¥µ}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡ÁÕ‰}å•…È¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}åÈ°å¤øÁÕ‰}å•…È°€‰%4¸àÁÕ‰±¥…Ñ¥½¸å•…Èˆ¤4(4(Œ€ä¸]•¥¡ÐÍ¡•µ”€¡Ý•¥¡Ñ•ÙÌÕ¹Ý•¥¡Ñ•¤4)¥µ}­}ÝÐ€ð´¥µ}¬€”ø”4(€™¥±Ñ•È¡ÝÑ}É½ÕÀ€•¥¸”Œ ‰Ý•¥¡Ñ•ˆ°€‰Õ¹Ý•¥¡Ñ•ˆ¤¤€”ø”4(€µÕÑ…Ñ”¡¥Í}Ý•¥¡Ñ•€ô¥™•±Í”¡ÝÑ}É½ÕÀ€ôô€‰Ý•¥¡Ñ•ˆ°€Ä°€À¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}ÝÐ°å¤ø¥Í}Ý•¥¡Ñ•°€‰%4¸äÝ•¥¡Ñ}Í¡•µ”€¡Ý•¥¡Ñ•ÙÌU\¤ˆ¤4(4(Œ€ÄÀ¸¹}É•…‘•ÉÌ€ ÈÙÌ€Ì¬ƒŠP™½È%4Ñ¡¥Ìµ…äÙ…ÉäÝ¥Ñ €ÄµÉ•…‘•È‘•Í¥¹Ì¤4)¥µ}­}¹È€ð´¥µ}¬€”ø”4(€™¥±Ñ•È …¥Ì¹¹„¡¹}É•…‘•ÉÌ¤¤€”ø”4(€µÕÑ…Ñ”¡¹É}…Ð€ô…Í•}Ý¡•¸ 4(€€€¹}É•…‘•ÉÌ€ôô€Äø€ˆÄˆ°4(€€€¹}É•…‘•ÉÌ€ôô€Èø€ˆÈˆ°4(€€€¹}É•…‘•ÉÌ€øô€Ìø€ˆÍÁ±ÕÌˆ°4(€€€QIUø9}¡…É…Ñ•É|4(€€¤¤€”ø”4(€™¥±Ñ•È …¥Ì¹¹„¡¹É}…Ð¤¤4)¥˜€¡±•¹Ñ ¡Õ¹¥ÅÕ”¡¥µ}­}¹È‘¹É}…Ð¤¤€øô€È¤ì4(€¥µ}­}¹È‘¹É}…Ð€ð´™…Ñ½È¡¥µ}­}¹È‘¹É}…Ð°±•Ù•±Ì€ôŒ ˆÄˆ°€ˆÈˆ°€ˆÍÁ±ÕÌˆ¤¤4(€ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}¹È°å¤ø¹É}…Ð°€‰%4¸ÄÀ¹}É•…‘•ÉÌ€ Ä¼È¼Ì¬¤ˆ¤4)ô4(4(Œ€ÄÄ¸EI0ÅÕ…±¥Ñä€¡1½ÜÙÌ5½‘•É…Ñ”½!¥ ¤4)¥µ}­}É½ˆ€ð´¥µ}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡¥Í}±½Ý}É½ˆ¤¤4)…Ð¡ÍÁÉ¥¹Ñ˜ ˆ€€€%4EI0ÅÕ…±¥Ñäè¬ô•€¡1½Üô•°5½‘•É…Ñ”½!¥ ô•¥q¸ˆ°4(€€€€€€€€€€€¹É½Ü¡¥µ}­}É½ˆ¤°ÍÕ´¡¥µ}­}É½ˆ‘¥Í}±½Ý}É½ˆ€ôô€Ä¤°4(€€€€€€€€€€€ÍÕ´¡¥µ}­}É½ˆ‘¥Í}±½Ý}É½ˆ€ôô€À¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}É½ˆ°å¤ø¥Í}±½Ý}É½ˆ°€‰%4¸ÄÄEI0ÅÕ…±¥Ñä€¡1½ÜÙÌ5½½!¥ ¤ˆ¤4(4(Œ€ÄÅˆ¸EI0å•Ìµ½Õ¹Ð€¡½¹Ñ¥¹Õ½ÕÌ¤4)¥µ}­}ÅÌ€ð´¥µ}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡Å…É•±}Í½É”¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}ÅÌ°å¤øÅ…É•±}Í½É”°€‰%4¸ÄÅˆEI0å•Ìµ½Õ¹Ð€¡½¹Ñ¥¹Õ½ÕÌ¤ˆ¤4(4(Œ€ÄÈ¸MÑÕ‘ä‘•Í¥¸€¡ÁÉ½ÍÁ•Ñ¥Ù”ÙÌÉ•ÑÉ½ÍÁ•Ñ¥Ù”¤4)¥µ}­}‘•Ì€ð´¥µ}¬€”ø”™¥±Ñ•È …¥Ì¹¹„¡¥Í}ÁÉ½ÍÁ•Ñ¥Ù”¤¤4)…Ð¡ÍÁÉ¥¹Ñ˜ ˆ€€€%4‘•Í¥¸‘…Ñ„…Ù…¥±…‰±”è¬ô•€¡É•ÑÉ¼ô•°ÁÉ½ÍÀô•¥q¸ˆ°4(€€€€€€€€€€€¹É½Ü¡¥µ}­}‘•Ì¤°ÍÕ´¡¥µ}­}‘•Ì‘¥Í}ÁÉ½ÍÁ•Ñ¥Ù”€ôô€À¤°4(€€€€€€€€€€€ÍÕ´¡¥µ}­}‘•Ì‘¥Í}ÁÉ½ÍÁ•Ñ¥Ù”€ôô€Ä¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}‘•Ì°å¤ø¥Í}ÁÉ½ÍÁ•Ñ¥Ù”°€‰%4¸ÄÈ‘•Í¥¸€¡ÁÉ½ÍÁ•Ñ¥Ù”ÙÌÉ•ÑÉ½ÍÁ•Ñ¥Ù”¤ˆ¤4(4(Œ€ÄÌ¸I•¥½¸€¡9½ÉÑ µ•É¥„€¼Í¥„€¼ÕÉ½Á”¤4)¥µ}­}É•œ€ð´¥µ}¬€”ø”4(€™¥±Ñ•È …¥Ì¹¹„¡É•¥½¹}…Ð¤¤€”ø”4(€µÕÑ…Ñ”¡É•¥½¹}…Ð€ô™…Ñ½È¡É•¥½¹}…Ð°±•Ù•±Ì€ôŒ ‰9½ÉÑ¡}µ•É¥„ˆ°€‰Í¥„ˆ°€‰ÕÉ½Á”ˆ¤¤¤4)…Ð¡ÍÁÉ¥¹Ñ˜ ˆ€€€%4É•¥½¸‘…Ñ„…Ù…¥±…‰±”è¬ô•€¡9ô•°Í¥„ô•°ÕÉ½Á”ô•¥q¸ˆ°4(€€€€€€€€€€€¹É½Ü¡¥µ}­}É•œ¤°4(€€€€€€€€€€€ÍÕ´¡¥µ}­}É•œ‘É•¥½¹}…Ð€ôô€‰9½ÉÑ¡}µ•É¥„ˆ¤°4(€€€€€€€€€€€ÍÕ´¡¥µ}­}É•œ‘É•¥½¹}…Ð€ôô€‰Í¥„ˆ¤°4(€€€€€€€€€€€ÍÕ´¡¥µ}­}É•œ‘É•¥½¹}…Ð€ôô€‰ÕÉ½Á”ˆ¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥µ}­}É•œ°å¤øÉ•¥½¹}…Ð°€‰%4¸ÄÌÉ•¥½¸€¡9½Í¥„½ÕÉ½Á”¤ˆ¤4(4(4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(Œ€Ô¸M9M%Q%Y%Qd€Äè]•¥¡Ñ•­…ÁÁ„½¹±ä€ ¬Pµ‘•É¥Ù•1\¤4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(4)…Ð ‰q¸ˆ°Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4)…Ð ‰M9M%Q%Y%Qd€Äè]•¥¡Ñ•­…ÁÁ„½¹±ä€ ¬Pµ‘•É¥Ù•1\¥q¸ˆ¤4)…Ð¡Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4(4(ŒAÉ•Á…É”Pµ‘•É¥Ù•1\‘…Ñ…Í•ÐÑ¡É½Õ Ñ¡”Í…µ”Á¥Á•±¥¹”…ÌÁÉ¥µ…Éä4)‘…Ñ}…Ñ}±Ü€ð´½µÀ€”ø”4(€™¥±Ñ•È¡•á±Õ‘”€ôô€‰…Ñ}‘•É¥Ù•‘}±Ý}Í•¹Í¥Ñ¥Ù¥Ñäˆ¤€”ø”4(€±•™Ñ}©½¥¸¡ÍÑÕ€”ø”Í•±•Ð¡ÍÑÕ‘å}¥°å•…È°ÁÕ‰±¥…Ñ¥½¹}ÑåÁ”°ÍÑÕ‘å}‘•Í¥¸¤°4(€€€€€€€€€€€‰ä€ô€‰ÍÑÕ‘å}¥ˆ¤€”ø”4(€±•™Ñ}©½¥¸¡Å…É•°€”ø”Í•±•Ð¡ÍÑÕ‘å}¥°½Ù•É…±±}É½ˆ°Å…É•±}å•Í}½Õ¹Ð¤°4(€€€€€€€€€€€‰ä€ô€‰ÍÑÕ‘å}¥ˆ¤€”ø”4(€µÕÑ…Ñ” 4(€€€­…ÁÁ„€ô…Ì¹¹Õµ•É¥Œ¡­…ÁÁ„¤°4(€€€¹Œ€ô…Ì¹¹Õµ•É¥Œ¡¹}åÍÑÌ¤°4(€€€ÍÑ…Ñ}ÑåÁ”€ô€‰­…ÁÁ„ˆ°4(€€€ÍÑ…Ð€ô­…ÁÁ„°4(€€€ÍÑ…Ñ}±…µÁ•€ôÁµ¥¸¡Áµ…à¡ÍÑ…Ð°€´À¸ääÔ¤°€À¸ääÔ¤°4(€€€å¤€ô…Ñ…¹ ¡ÍÑ…Ñ}±…µÁ•¤°4(€€€Ù¥}¸€ô¥™•±Í” …¥Ì¹¹„¡¹Œ¤€˜¹Œ€ø€Ì°€Ä€¼€¡¹Œ€´€Ì¤°9}É•…±|¤°4(€€€Ù¤€ôÙ¥}¸°4(€€€ÝÑ}É½ÕÀ€ô€‰Ý•¥¡Ñ•ˆ4(€€¤€”ø”4(€™¥±Ñ•È …¥Ì¹¹„¡å¤¤€˜€…¥Ì¹¹„¡Ù¤¤€˜Ù¤€ø€À¤4)…Ð¡ÍÁÉ¥¹Ñ˜ ‰Pµ‘•É¥Ù•1\É½ÝÌ±½…‘•è€•‘q¸ˆ°¹É½Ü¡‘…Ñ}…Ñ}±Ü¤¤¤4(4)…Ñ}±Ý}¥È€ð´‘…Ñ}…Ñ}±Ü€”ø”™¥±Ñ•È¡…¹…±åÑ¥}ÍÑÉ…ÑÕ´€ôô€‰¥¹Ñ•ÈµÉ•…‘•Èˆ¤4)…Ñ}±Ý}¥´€ð´‘…Ñ}…Ñ}±Ü€”ø”™¥±Ñ•È¡…¹…±åÑ¥}ÍÑÉ…ÑÕ´€ôô€‰¥¹Ñ•Èµµ½‘…±¥Ñäˆ¤4(4(Œ]•¥¡Ñ•Í•¹Í¥Ñ¥Ù¥Ñäè…ÕÑ¡½ÈµÉ•Á½ÉÑ•Ý•¥¡Ñ•€¬Pµ‘•É¥Ù•1\4)¥É}Ý¬€ð´‰¥¹‘}É½ÝÌ 4(€¥É}¬€”ø”™¥±Ñ•È¡ÝÑ}É½ÕÀ€ôô€‰Ý•¥¡Ñ•ˆ¤°4(€…Ñ}±Ý}¥È4(¤4)¥É}Ý­}É•Ì€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥É}Ý¬°€‰%HÝ•¥¡Ñ•­…ÁÁ„€ ­Pµ‘•É¥Ù•1\¤ˆ¤4)…‘‘}É½Ü ‰%HÝ•¥¡Ñ•IYˆ°¥É}Ý­}É•Ì°¥É}Ý­}É•Ì‘¬°¥É}Ý­}É•Ì‘¹}±ÕÍÑ•ÉÌ¤4(4)¥µ}Ý¬€ð´‰¥¹‘}É½ÝÌ 4(€¥µ}¬€”ø”™¥±Ñ•È¡ÝÑ}É½ÕÀ€ôô€‰Ý•¥¡Ñ•ˆ¤°4(€…Ñ}±Ý}¥´4(¤4)¥µ}Ý­}É•Ì€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥µ}Ý¬°€‰%4Ý•¥¡Ñ•­…ÁÁ„€ ­Pµ‘•É¥Ù•1\¤ˆ¤4(4(4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(Œ€Ø¸M9M%Q%Y%Qd€ÈèU¹Ý•¥¡Ñ•­…ÁÁ„½¹±ä4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(4)…Ð ‰q¸ˆ°Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4)…Ð ‰M9M%Q%Y%Qd€ÈèU¹Ý•¥¡Ñ•­…ÁÁ„½¹±åq¸ˆ¤4)…Ð¡Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4(4)¥É}Õ¬€ð´¥É}¬€”ø”™¥±Ñ•È¡ÝÑ}É½ÕÀ€ôô€‰Õ¹Ý•¥¡Ñ•ˆ¤4)¥É}Õ­}É•Ì€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥É}Õ¬°€‰%HÕ¹Ý•¥¡Ñ•­…ÁÁ„ˆ¤4)…‘‘}É½Ü ‰%HÕ¹Ý•¥¡Ñ•IYˆ°¥É}Õ­}É•Ì°¥É}Õ­}É•Ì‘¬°¥É}Õ­}É•Ì‘¹}±ÕÍÑ•ÉÌ¤4(4)¥µ}Õ¬€ð´¥µ}¬€”ø”™¥±Ñ•È¡ÝÑ}É½ÕÀ€ôô€‰Õ¹Ý•¥¡Ñ•ˆ¤4)¥µ}Õ­}É•Ì€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥µ}Õ¬°€‰%4Õ¹Ý•¥¡Ñ•­…ÁÁ„ˆ¤4)…‘‘}É½Ü ‰%4Õ¹Ý•¥¡Ñ•IYˆ°¥µ}Õ­}É•Ì°¥µ}Õ­}É•Ì‘¬°¥µ}Õ­}É•Ì‘¹}±ÕÍÑ•ÉÌ¤4(4(4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(Œ€Ü¸M9M%Q%Y%Qd€Ìè€Ôµ…Ñ•½Éä­…ÁÁ„½¹±ä4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(4)…Ð ‰q¸ˆ°Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4)…Ð ‰M9M%Q%Y%Qd€Ìè€Ôµ…Ñ•½Éä­…ÁÁ„½¹±åq¸ˆ¤4)…Ð¡Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4(4)¥É|ÕŒ€ð´¥É}¬€”ø”™¥±Ñ•È¡¹}…Ñ•½É¥•Ì€ôô€Ô¤4)¥É|Õ}É•Ì€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥É|ÕŒ°€‰%H€Ôµ…Ð­…ÁÁ„ˆ¤4)…‘‘}É½Ü ‰%H€Ôµ…ÐIYˆ°¥É|Õ}É•Ì°¥É|Õ}É•Ì‘¬°¥É|Õ}É•Ì‘¹}±ÕÍÑ•ÉÌ¤4(4)¥µ|ÕŒ€ð´¥µ}¬€”ø”™¥±Ñ•È¡¹}…Ñ•½É¥•Ì€ôô€Ô¤4)¥µ|Õ}É•Ì€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥µ|ÕŒ°€‰%4€Ôµ…Ð­…ÁÁ„ˆ¤4)…‘‘}É½Ü ‰%4€Ôµ…ÐIYˆ°¥µ|Õ}É•Ì°¥µ|Õ}É•Ì‘¬°¥µ|Õ}É•Ì‘¹}±ÕÍÑ•ÉÌ¤4(4(4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(Œ€à¸aA1=IQ=Idè±°µ•ÑÉ¥ÌÁ½½±•4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(4)…Ð ‰q¸ˆ°Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4)…Ð ‰aA1=IQ=Idè±°µ•ÑÉ¥ÌÁ½½±•€¡­…ÁÁ„€¬Ý•Ð€¬%¥q¸ˆ¤4)…Ð¡Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4(4)¥É}…±°€ð´‘…Ð€”ø”™¥±Ñ•È¡…¹…±åÑ¥}ÍÑÉ…ÑÕ´€ôô€‰¥¹Ñ•ÈµÉ•…‘•Èˆ¤4)¥É}…±±}É•Ì€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥É}…±°°€‰%H…±°µ•ÑÉ¥Ìˆ¤4)…‘‘}É½Ü ‰%H…±°µµ•ÑÉ¥ŒIYˆ°¥É}…±±}É•Ì°¥É}…±±}É•Ì‘¬°¥É}…±±}É•Ì‘¹}±ÕÍÑ•ÉÌ¤4(4(ŒÍÑ…Ñ}ÑåÁ”µ½‘•É…Ñ½È4)¥É}…±±}ÍÐ€ð´¥É}…±°€”ø”4(€µÕÑ…Ñ”¡ÍÑ…Ñ}ÑåÁ”€ô™…Ñ½È¡ÍÑ…Ñ}ÑåÁ”°±•Ù•±Ì€ôŒ ‰­…ÁÁ„ˆ°€‰Ý•Ðˆ°€‰¥Œˆ¤¤¤4)ÉÕ¹}µ½‘•É…Ñ½È¡¥É}…±±}ÍÐ°å¤øÍÑ…Ñ}ÑåÁ”°€‰ÍÑ…Ñ}ÑåÁ”€¡­…ÁÁ„½Ý•Ð½¥Œ¤ˆ¤4(4)¥µ}…±°€ð´‘…Ð€”ø”™¥±Ñ•È¡…¹…±åÑ¥}ÍÑÉ…ÑÕ´€ôô€‰¥¹Ñ•Èµµ½‘…±¥Ñäˆ¤4)¥µ}…±±}É•Ì€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥µ}…±°°€‰%4…±°µ•ÑÉ¥Ìˆ¤4)…‘‘}É½Ü ‰%4…±°µµ•ÑÉ¥ŒIYˆ°¥µ}…±±}É•Ì°¥µ}…±±}É•Ì‘¬°¥µ}…±±}É•Ì‘¹}±ÕÍÑ•ÉÌ¤4(4(4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(Œ€ä¸aA1=IQ=IdèMµ…±°ÍÑÉ…Ñ„…±°µ•™™•ÑÌ€¡­…ÁÁ„µ½¹±ä¤4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(4)…Ð ‰q¸ˆ°Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4)…Ð ‰aA1=IQ=IdèMµ…±°ÍÑÉ…Ñ„…±°µ•™™•ÑÌ€¡­…ÁÁ„µ½¹±ä¥q¸ˆ¤4)…Ð¡Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4(4)¥¹ÑÉ…}…±±}¬€ð´‘…Ñ}­…ÁÁ„€”ø”™¥±Ñ•È¡…¹…±åÑ¥}ÍÑÉ…ÑÕ´€ôô€‰¥¹ÑÉ„µÉ•…‘•Èˆ¤4)¥¹ÑÉ…}ÉÙ”€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥¹ÑÉ…}…±±}¬°€‰%¹ÑÉ„µÉ•…‘•È­…ÁÁ„…±°µ•™™•ÑÌˆ¤4(4)¥Ù}…±±}¬€ð´‘…Ñ}­…ÁÁ„€”ø”™¥±Ñ•È¡…¹…±åÑ¥}ÍÑÉ…ÑÕ´€ôô€‰¥¹Ñ•ÈµÙ•ÉÍ¥½¸ˆ¤4)¥Ù}ÉÙ”€ð´ÉÕ¹}ÉÙ•}ÈÈ¡¥Ù}…±±}¬°€‰%¹Ñ•ÈµÙ•ÉÍ¥½¸­…ÁÁ„…±°µ•™™•ÑÌˆ¤4(4(4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(Œ€ÄÀ¸IYÙÌ=AL=5AI%M=84(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(4)…Ð ‰q¸ˆ°Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4)…Ð ‰IYÙÌ=AL=5AI%M=8€¡­…ÁÁ„µ½¹±ä¥q¸ˆ¤4)…Ð¡Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4(4)½µÁ…É”€ð´™Õ¹Ñ¥½¸¡ÉÙ”°½Á•Ì°±…‰•°¤ì4(€¥˜€¡¥Ì¹¹Õ±°¡ÉÙ”¤ñð¥Ì¹¹Õ±°¡½Á•Ì¤¤É•ÑÕÉ¸ ¤4(€‘¥™˜€ð´€¡½Á•Ì‘Á½½±•‘}­…ÁÁ„€´ÉÙ”‘Á½½±•‘}­…ÁÁ„¤€¼ÉÙ”‘Á½½±•‘}­…ÁÁ„€¨€ÄÀÀ4(€½Ù•É±…À€ð´€„¡½Á•Ì‘¥}­…ÁÁ…lÉt€ðÉÙ”‘¥}­…ÁÁ…lÅtðÉÙ”‘¥}­…ÁÁ…lÉt€ð½Á•Ì‘¥}­…ÁÁ…lÅt¤4(€…Ð¡ÍÁÉ¥¹Ñ˜ ‰q¸•Ìéq¸€IYè€”¸Í˜€ ”¸Í˜°€”¸Í˜¥q¸€=ALè€”¸Í˜€ ”¸Í˜°€”¸Í˜¥q¸€¥™˜è€”¬¸Å˜””°$½Ù•É±…Àè€•Íq¸ˆ°4(€€€€€€€€€€€€€±…‰•°°4(€€€€€€€€€€€€€ÉÙ”‘Á½½±•‘}­…ÁÁ„°ÉÙ”‘¥}­…ÁÁ…lÅt°ÉÙ”‘¥}­…ÁÁ…lÉt°4(€€€€€€€€€€€€€½Á•Ì‘Á½½±•‘}­…ÁÁ„°½Á•Ì‘¥}­…ÁÁ…lÅt°½Á•Ì‘¥}­…ÁÁ…lÉt°4(€€€€€€€€€€€€€‘¥™˜°¥™•±Í”¡½Ù•É±…À°€‰eLˆ°€‰9<ˆ¤¤¤4)ô4(4)½µÁ…É”¡¥É}É•Ì°¥É}½Á•Í}É•Ì°€‰%¹Ñ•ÈµÉ•…‘•Èˆ¤4)½µÁ…É”¡¥µ}É•Ì°¥µ}½Á•Í}É•Ì°€‰%¹Ñ•Èµµ½‘…±¥Ñäˆ¤4(4(4(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(Œ€ÄÄ¸MU55IdQ	14(Œ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(4)…Ð ‰q¸ˆ°Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4)…Ð ‰MU55IdQ	1q¸ˆ¤4)…Ð¡Á…ÍÑ”¡É•À ˆôˆ°€ÜÀ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4(4)…Ð¡ÍÁÉ¥¹Ñ˜ ‰q¸”´ÈÕÌ€”ÕÌ€”ÙÌ€”áÌ€”ÄáÌ€”ÄáÍq¸ˆ°4(€€€€€€€€€€€€‰¹…±åÍ¥Ìˆ°€‰¬ˆ°€‰±ÕÍÐˆ°€‰•ÍÐˆ°€ˆäÔ”$ˆ°€ˆäÔ”A$ˆ¤¤4)…Ð¡Á…ÍÑ”¡É•À ˆ´ˆ°€àÔ¤°½±±…ÁÍ”€ô€ˆˆ¤°€‰q¸ˆ¤4(4)™½È€¡È¥¸ÍÕµµ…Éå}É½ÝÌ¤ì4(€…Ð¡ÍÁÉ¥¹Ñ˜ ˆ”´ÈÕÌ€”Õ€”Ù€”à¸Í˜€l”¸Í˜°€”¸Í™t€l”¸Í˜°€”¸Í™uq¸ˆ°4(€€€€€€€€€€€€€È‘±…‰•°°È‘¬°È‘±ÕÍÐ°È‘­…ÁÁ„°4(€€€€€€€€€€€€€È‘¥lÅt°È‘¥lÉt°È‘Á¥lÅt°È‘Á¥lÉt¤¤4)ô4(4)…Ð ‰q¹½¹”¹q¸ˆ¤4(