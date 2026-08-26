## Reproduce the descriptive inter-modality management-tier synthesis from the
## public per-comparison 2x2 reconstructions.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 0L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
repo_dir <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/")
out_dir <- file.path(repo_dir, "revision", "analysis", "results")
input_path <- file.path(out_dir, "intermodality_management_tier_rows.csv")

if (!file.exists(input_path)) stop("Missing input: ", input_path)
rows <- read.csv(
  input_path,
  check.names = FALSE,
  na.strings = c("", "NA"),
  fileEncoding = "UTF-8-BOM"
)

required <- c(
  "tier_row_id", "modality_pair", "study", "n", "lower_lower",
  "other_higher", "other_lower", "higher_higher", "observed_agreement",
  "kappa", "ac1", "net_other_higher", "reconstruction_status"
)
missing_columns <- setdiff(required, names(rows))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

expected_ids <- sprintf("IMT%02d", 1:26)
stopifnot(
  nrow(rows) == 26L,
  identical(rows$tier_row_id, expected_ids),
  !anyDuplicated(rows$tier_row_id),
  all(rows$modality_pair %in% c("CT-MRI", "CT-CEUS")),
  all(rows$reconstruction_status %in% c("full", "direction-only"))
)

full <- rows[rows$reconstruction_status == "full", ]
direction_only <- rows[rows$reconstruction_status == "direction-only", ]
stopifnot(nrow(full) == 25L, nrow(direction_only) == 1L)

cell_columns <- c("lower_lower", "other_higher", "other_lower", "higher_higher")
stopifnot(
  all(complete.cases(full[, cell_columns])),
  all(complete.cases(full[, c("observed_agreement", "kappa", "ac1")])),
  all(full$net_other_higher == full$other_higher - full$other_lower)
)

recompute_measures <- function(x) {
  n <- x$lower_lower + x$other_higher + x$other_lower + x$higher_higher
  observed <- (x$lower_lower + x$higher_higher) / n
  ct_lower <- (x$lower_lower + x$other_higher) / n
  ct_higher <- 1 - ct_lower
  other_lower <- (x$lower_lower + x$other_lower) / n
  other_higher <- 1 - other_lower
  expected_kappa <- ct_lower * other_lower + ct_higher * other_higher
  kappa <- (observed - expected_kappa) / (1 - expected_kappa)
  mean_lower <- (ct_lower + other_lower) / 2
  expected_ac1 <- 2 * mean_lower * (1 - mean_lower)
  ac1 <- (observed - expected_ac1) / (1 - expected_ac1)
  c(observed_agreement = observed, kappa = kappa, ac1 = ac1)
}

recomputed <- t(vapply(
  seq_len(nrow(full)),
  function(i) recompute_measures(full[i, , drop = FALSE]),
  numeric(3)
))
stopifnot(
  max(abs(recomputed[, "observed_agreement"] - full$observed_agreement)) < 5e-4,
  max(abs(recomputed[, "kappa"] - full$kappa)) < 5e-4,
  max(abs(recomputed[, "ac1"] - full$ac1)) < 5e-4
)

summarize_rows <- function(x, status) {
  data.frame(
    modality_pair = unique(x$modality_pair),
    reconstruction_status = status,
    comparisons = nrow(x),
    studies = length(unique(x$study)),
    lesions_total = sum(x$n),
    observed_agreement_min = min(x$observed_agreement),
    observed_agreement_max = max(x$observed_agreement),
    kappa_min = if (all(is.na(x$kappa))) NA_real_ else min(x$kappa, na.rm = TRUE),
    kappa_max = if (all(is.na(x$kappa))) NA_real_ else max(x$kappa, na.rm = TRUE),
    ac1_min = if (all(is.na(x$ac1))) NA_real_ else min(x$ac1, na.rm = TRUE),
    ac1_max = if (all(is.na(x$ac1))) NA_real_ else max(x$ac1, na.rm = TRUE),
    net_upward_comparisons = sum(x$net_other_higher > 0),
    net_balanced_comparisons = sum(x$net_other_higher == 0),
    net_downward_comparisons = sum(x$net_other_higher < 0),
    upward_crossings = sum(x$other_higher),
    downward_crossings = sum(x$other_lower)
  )
}

summary_rows <- rbind(
  summarize_rows(full[full$modality_pair == "CT-MRI", ], "full"),
  summarize_rows(full[full$modality_pair == "CT-CEUS", ], "full"),
  summarize_rows(direction_only, "direction-only")
)

stopifnot(
  identical(summary_rows$comparisons, c(15L, 10L, 1L)),
  identical(summary_rows$studies, c(10L, 9L, 1L)),
  identical(summary_rows$net_upward_comparisons, c(13L, 7L, 0L)),
  identical(summary_rows$net_balanced_comparisons, c(1L, 1L, 0L)),
  identical(summary_rows$net_downward_comparisons, c(1L, 2L, 1L))
)

output_path <- file.path(out_dir, "intermodality_management_tier_summary.csv")
write.csv(summary_rows, output_path, row.names = FALSE, na = "")

print(summary_rows)
message("Wrote ", output_path)
