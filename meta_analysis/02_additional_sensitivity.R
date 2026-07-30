## ============================================================
## Bosniak SRMA — Additional Sensitivity Analyses
##
## 1. Subgroup pooled estimates (OPES, by modality / modality pair)
## 2. Direct-reported only sensitivity (excluding synthetic/restored/derived)
## 3. Classified-weight-only sensitivity (excluding weight_unspecified)
## ============================================================

library(metafor)
library(clubSandwich)
library(dplyr)
library(readr)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 0) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
setwd(script_dir)

comp <- read_csv("comparisons.csv", show_col_types = FALSE)
stud <- read_csv("studies.csv", show_col_types = FALSE)
comp <- comp %>%
  left_join(stud %>% select(study_id, year, publication_type), by = "study_id")

# --- Data prep (same as 01_main_analysis.R) ---
dat <- comp %>%
  filter(is.na(exclude) | exclude == "") %>%
  mutate(
    kappa = as.numeric(kappa),
    stat = kappa,
    ci_lo = as.numeric(ci_lower), ci_hi = as.numeric(ci_upper),
    se_raw = as.numeric(kappa_se),
    nc = as.numeric(n_cysts),
    stat_clamped = pmin(pmax(stat, -0.995), 0.995),
    yi = atanh(stat_clamped),
    vi_ci = ifelse(!is.na(ci_lo) & !is.na(ci_hi),
      ((atanh(pmin(pmax(ci_hi,-0.995),0.995)) - atanh(pmin(pmax(ci_lo,-0.995),0.995))) /
        (2*qnorm(0.975)))^2, NA_real_),
    vi_se = ifelse(!is.na(se_raw), (se_raw/(1-stat_clamped^2))^2, NA_real_),
    vi_n = ifelse(!is.na(nc) & nc > 3, 1/(nc-3), NA_real_),
    vi = coalesce(vi_ci, vi_se, vi_n),
    wt_group = case_when(
      weight_scheme %in% c("linear","quadratic","custom","Cicchetti","weighted_unspecified") ~ "weighted",
      weight_scheme == "unweighted" ~ "unweighted",
      TRUE ~ "other"
    ),
    mod_group = case_when(
      std_modality_1 == "CT" ~ "CT",
      std_modality_1 == "MRI" ~ "MRI",
      std_modality_1 %in% c("CEUS","US","SMI") ~ "CEUS_US",
      TRUE ~ "other"
    ),
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
  filter(!is.na(stat) & !is.na(vi) & vi > 0)

dat_kappa <- dat %>% filter(!is.na(kappa))


# ============================================================
# Helper: RVE+CR2
# ============================================================
run_rve <- function(data, label) {
  if (nrow(data) < 3 || length(unique(data$dep_id)) < 3) {
    cat(sprintf("  %s: SKIPPED (k=%d, clusters=%d)\n", label, nrow(data), length(unique(data$dep_id))))
    return(NULL)
  }
  mod <- rma.mv(yi, vi, random = list(~1|dep_id, ~1|comparison_id), data = data, method = "REML")
  cr2 <- coef_test(mod, vcov = "CR2", cluster = data$dep_id)
  tc <- qt(0.975, cr2$df_Satt)
  k <- tanh(mod$b[1])
  ci <- tanh(c(mod$b[1] - tc*cr2$SE, mod$b[1] + tc*cr2$SE))
  cat(sprintf("  %s: k=%d, clusters=%d, kappa=%.3f (%.3f, %.3f), df=%.1f\n",
              label, nrow(data), length(unique(data$dep_id)), k, ci[1], ci[2], cr2$df_Satt))
  list(kappa = k, ci = ci, k_eff = nrow(data), clusters = length(unique(data$dep_id)))
}

# Helper: OPES RE
run_opes <- function(data, label) {
  if (nrow(data) < 2) {
    cat(sprintf("  %s: SKIPPED (k=%d)\n", label, nrow(data)))
    return(NULL)
  }
  data <- data %>%
    mutate(
      opes_val_str = ifelse(
        grepl("^(median_synthetic|reported_overall_restored):", opes_basis),
        sub("^[^:]+:", "", opes_basis), NA_character_),
      opes_val = coalesce(as.numeric(opes_val_str), stat),
      opes_clamped = pmin(pmax(opes_val, -0.995), 0.995),
      yi_o = atanh(opes_clamped),
      vi_o = case_when(
        grepl("^(median_synthetic|reported_overall_restored):", opes_basis) ~
          ifelse(!is.na(nc) & nc > 3, 1/(nc-3), vi),
        TRUE ~ vi
      )
    )
  mod <- rma(yi_o, vi_o, data = data, method = "REML")
  k <- tanh(mod$b[1])
  ci <- tanh(c(mod$ci.lb, mod$ci.ub))
  cat(sprintf("  %s: k=%d, kappa=%.3f (%.3f, %.3f), I2=%.1f%%\n",
              label, nrow(data), k, ci[1], ci[2], mod$I2))
  list(kappa = k, ci = ci, k_eff = nrow(data), I2 = mod$I2)
}


# ============================================================
# 1. SUBGROUP POOLED ESTIMATES (OPES)
# ============================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("1. SUBGROUP POOLED ESTIMATES (OPES)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

opes <- dat_kappa %>% filter(opes_include == 1)

# IR by modality
cat("\n--- Inter-reader by modality ---\n")
ir_opes <- opes %>% filter(analytic_stratum == "inter-reader")
for (g in c("CT", "MRI", "CEUS_US", "other")) {
  sub <- ir_opes %>% filter(mod_group == g)
  run_opes(sub, paste0("IR ", g))
}
cat("  IR overall:\n")
run_opes(ir_opes, "IR all")

# IM by modality pair
cat("\n--- Inter-modality by modality pair ---\n")
im_opes <- opes %>% filter(analytic_stratum == "inter-modality")
for (g in c("CT_MRI", "CT_CEUS", "MRI_CEUS", "other")) {
  sub <- im_opes %>% filter(mod_pair == g)
  run_opes(sub, paste0("IM ", g))
}
cat("  IM overall:\n")
run_opes(im_opes, "IM all")


# ============================================================
# 2. DIRECT-REPORTED ONLY SENSITIVITY
# ============================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("2. DIRECT-REPORTED ONLY SENSITIVITY\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\n--- RVE+CR2 (all kappa effects, excluding derived_from_table) ---\n")
ir_direct <- dat_kappa %>%
  filter(analytic_stratum == "inter-reader",
         is.na(metric_note) | !grepl("derived_from_table", metric_note))
run_rve(ir_direct, "IR direct-reported RVE")

im_direct <- dat_kappa %>%
  filter(analytic_stratum == "inter-modality",
         is.na(metric_note) | !grepl("derived_from_table", metric_note))
run_rve(im_direct, "IM direct-reported RVE")

cat("\n--- OPES (excluding synthetic/restored/derived) ---\n")
ir_opes_direct <- ir_opes %>%
  filter(!grepl("^(median_synthetic|reported_overall_restored|derived_from_table)", opes_basis),
         is.na(metric_note) | !grepl("derived_from_table", metric_note))
run_opes(ir_opes_direct, "IR OPES direct-only")

im_opes_direct <- im_opes %>%
  filter(!grepl("^(median_synthetic|reported_overall_restored|derived_from_table)", opes_basis),
         is.na(metric_note) | !grepl("derived_from_table", metric_note))
run_opes(im_opes_direct, "IM OPES direct-only")

# Show what was excluded
cat("\n  Excluded from IR OPES:\n")
ir_excl <- ir_opes %>% filter(grepl("^(median_synthetic|reported_overall_restored)", opes_basis))
if (nrow(ir_excl) > 0) {
  for (i in 1:nrow(ir_excl)) cat(sprintf("    %s (%s): basis=%s\n", ir_excl$comparison_id[i], ir_excl$study_id[i], ir_excl$opes_basis[i]))
} else cat("    none\n")

cat("  Excluded from IM OPES:\n")
im_excl <- im_opes %>% filter(grepl("^(median_synthetic|reported_overall_restored)", opes_basis))
if (nrow(im_excl) > 0) {
  for (i in 1:nrow(im_excl)) cat(sprintf("    %s (%s): basis=%s\n", im_excl$comparison_id[i], im_excl$study_id[i], im_excl$opes_basis[i]))
} else cat("    none\n")

# Derived from table in all effects
cat("\n  derived_from_table effects:\n")
dft <- dat_kappa %>% filter(!is.na(metric_note), grepl("derived_from_table", metric_note))
cat(sprintf("    Total: %d (IR: %d, IM: %d)\n",
    nrow(dft),
    sum(dft$analytic_stratum == "inter-reader"),
    sum(dft$analytic_stratum == "inter-modality")))


# ============================================================
# 3. CLASSIFIED-WEIGHT-ONLY SENSITIVITY
# ============================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("3. CLASSIFIED-WEIGHT-ONLY SENSITIVITY (excl. unspecified)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\n--- RVE+CR2 (weighted + unweighted, excluding unspecified) ---\n")
ir_cw <- dat_kappa %>%
  filter(analytic_stratum == "inter-reader",
         wt_group %in% c("weighted", "unweighted"))
run_rve(ir_cw, "IR classified-weight RVE")

im_cw <- dat_kappa %>%
  filter(analytic_stratum == "inter-modality",
         wt_group %in% c("weighted", "unweighted"))
run_rve(im_cw, "IM classified-weight RVE")

cat("\n--- OPES (excluding weight_unspecified) ---\n")
ir_opes_cw <- ir_opes %>%
  filter(weight_scheme != "weighted_unspecified" | is.na(weight_scheme))
run_opes(ir_opes_cw, "IR OPES classified-wt")

im_opes_cw <- im_opes %>%
  filter(weight_scheme != "weighted_unspecified" | is.na(weight_scheme))
run_opes(im_opes_cw, "IM OPES classified-wt")

# Count excluded
cat(sprintf("\n  Excluded unspecified: IR RVE %d, IM RVE %d, IR OPES %d, IM OPES %d\n",
    nrow(dat_kappa %>% filter(analytic_stratum=="inter-reader", wt_group=="other")),
    nrow(dat_kappa %>% filter(analytic_stratum=="inter-modality", wt_group=="other")),
    sum(ir_opes$weight_scheme == "weighted_unspecified", na.rm=TRUE),
    sum(im_opes$weight_scheme == "weighted_unspecified", na.rm=TRUE)))


cat("\nDone.\n")
