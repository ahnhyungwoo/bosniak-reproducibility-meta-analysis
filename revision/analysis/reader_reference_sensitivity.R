## Sensitivity: add back reader-versus-reference effects.
##
## The primary inter-reader set holds comparisons between readers who
## classified the same lesions independently. Effects that compare a reader
## with a consensus or gold-standard reading carry
## analytic_stratum == "reader-reference" in comparisons.csv and are kept out
## of that set. This script reports the primary model, the model with those
## effects added back (the set analysed in the first revision), and a
## description of the reader-versus-reference effects themselves.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
project_dir <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/")
local_lib <- file.path(project_dir, "revision", "R_lib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(metafor)
  library(clubSandwich)
})

comp <- read_csv(file.path(project_dir, "meta_analysis", "comparisons.csv"), show_col_types = FALSE)

dat <- comp %>%
  filter(is.na(exclude) | trimws(exclude) == "") %>%
  mutate(
    kappa_num = suppressWarnings(as.numeric(kappa)),
    ci_lo = suppressWarnings(as.numeric(ci_lower)),
    ci_hi = suppressWarnings(as.numeric(ci_upper)),
    se_raw = suppressWarnings(as.numeric(kappa_se)),
    nc = suppressWarnings(as.numeric(n_cysts)),
    kappa_clamped = pmin(pmax(kappa_num, -0.995), 0.995),
    yi_z = atanh(kappa_clamped),
    vi_z_ci = if_else(
      !is.na(ci_lo) & !is.na(ci_hi),
      ((atanh(pmin(pmax(ci_hi, -0.995), 0.995)) -
          atanh(pmin(pmax(ci_lo, -0.995), 0.995))) / (2 * qnorm(0.975)))^2,
      NA_real_
    ),
    vi_z_se = if_else(!is.na(se_raw), (se_raw / (1 - kappa_clamped^2))^2, NA_real_),
    vi_z_n = if_else(!is.na(nc) & nc > 3, 1 / (nc - 3), NA_real_),
    vi_z = coalesce(vi_z_ci, vi_z_se, vi_z_n)
  ) %>%
  filter(
    analytic_stratum %in% c("inter-reader", "reader-reference"),
    !is.na(kappa_num), !is.na(vi_z), vi_z > 0
  )

direct <- dat %>% filter(analytic_stratum == "inter-reader")
reference <- dat %>% filter(analytic_stratum == "reader-reference")
if (nrow(reference) == 0) stop("No reader-versus-reference effects found.")
if (any(reference$reader_structure != "vs_reference")) {
  stop("Reader-versus-reference effects must carry reader_structure vs_reference.")
}

run_rve <- function(d, label) {
  model <- rma.mv(
    yi_z, vi_z,
    random = list(~ 1 | dep_id, ~ 1 | comparison_id),
    data = d,
    method = "REML"
  )
  robust <- coef_test(model, vcov = "CR2", cluster = d$dep_id)
  estimate_z <- robust$beta[1]
  se_z <- robust$SE[1]
  df <- robust$df_Satt[1]
  critical <- qt(0.975, df)
  tibble(
    analysis = label,
    effects = nrow(d),
    studies = n_distinct(d$study_id),
    clusters = n_distinct(d$dep_id),
    estimate = tanh(estimate_z),
    ci_low = tanh(estimate_z - critical * se_z),
    ci_high = tanh(estimate_z + critical * se_z),
    df = df,
    p = robust$p_Satt[1],
    pi_low = tanh(estimate_z - critical * sqrt(se_z^2 + sum(model$sigma2))),
    pi_high = tanh(estimate_z + critical * sqrt(se_z^2 + sum(model$sigma2)))
  )
}

results <- bind_rows(
  run_rve(direct, "Direct reader comparisons (primary)"),
  run_rve(dat, "Including reader-versus-reference effects")
)

reference_by_study <- reference %>%
  group_by(study_id) %>%
  summarise(
    effects = n(),
    comparison_ids = paste(comparison_id, collapse = "; "),
    kappa_min = min(kappa_num),
    kappa_median = median(kappa_num),
    kappa_max = max(kappa_num),
    .groups = "drop"
  )

out_dir <- file.path(project_dir, "revision", "analysis", "results")
write_csv(results, file.path(out_dir, "reader_reference_sensitivity.csv"))
write_csv(reference_by_study, file.path(out_dir, "reader_reference_effects_by_study.csv"))
print(results, width = Inf)
print(reference_by_study, width = Inf)
