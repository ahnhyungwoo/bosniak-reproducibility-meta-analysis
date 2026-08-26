suppressPackageStartupMessages({
  library(metafor)
  library(clubSandwich)
})

options(stringsAsFactors = FALSE)

out_dir <- file.path("revision", "analysis", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

records <- list()

add_cells <- function(study_id, comparison, source_type, ll, lh, hl, hh,
                      scale = "Bosniak", notes = "") {
  n <- ll + lh + hl + hh
  records[[length(records) + 1L]] <<- data.frame(
    study_id = study_id,
    comparison = comparison,
    source_type = source_type,
    scale = scale,
    lower_lower = ll,
    lower_higher = lh,
    higher_lower = hl,
    higher_higher = hh,
    n = n,
    included = 1L,
    exclusion_reason = "",
    notes = notes
  )
}

add_excluded <- function(study_id, comparison, source_type, n, scale, reason, notes = "") {
  records[[length(records) + 1L]] <<- data.frame(
    study_id = study_id,
    comparison = comparison,
    source_type = source_type,
    scale = scale,
    lower_lower = NA_integer_,
    lower_higher = NA_integer_,
    higher_lower = NA_integer_,
    higher_higher = NA_integer_,
    n = n,
    included = 0L,
    exclusion_reason = reason,
    notes = notes
  )
}

collapse_matrix <- function(x, low_rows, low_cols = low_rows) {
  high_rows <- setdiff(seq_len(nrow(x)), low_rows)
  high_cols <- setdiff(seq_len(ncol(x)), low_cols)
  c(
    ll = sum(x[low_rows, low_cols, drop = FALSE]),
    lh = sum(x[low_rows, high_cols, drop = FALSE]),
    hl = sum(x[high_rows, low_cols, drop = FALSE]),
    hh = sum(x[high_rows, high_cols, drop = FALSE])
  )
}

add_matrix <- function(study_id, comparison, source_type, x, low_rows, low_cols = low_rows,
                       scale = "Bosniak", notes = "") {
  cc <- collapse_matrix(x, low_rows, low_cols)
  add_cells(study_id, comparison, source_type, cc["ll"], cc["lh"], cc["hl"], cc["hh"],
            scale, notes)
}

# Directly reported reader-reader or reader-reader-like cross-tabulations.
add_cells("Benjaminov_2006", "Reviewer 1 vs Reviewer 2", "direct_reader_pair", 5, 2, 2, 23)
add_cells("Bertolotto_2015", "Radiologist 2 vs Radiologist 1", "direct_reader_pair", 11, 1, 4, 7,
          notes = "Study sample contained categories IIF-IV.")
add_cells("El-Mokadem_2014", "Initial report vs multidisciplinary reassessment", "clinical_reassessment", 56, 5, 7, 28)
add_cells("Kim_2010", "Reader 2 vs Reader 1", "direct_reader_pair", 62, 4, 13, 46)
add_cells("Lucocq_2021", "Initial report vs specialist reassessment", "clinical_reassessment", 153, 18, 15, 71,
          notes = "One raw cross-tabulation underlies two weighted-kappa comparison rows.")
add_cells("Lucocq_2024", "Reader 1 vs Reader 2, version 2005", "direct_reader_pair", 6, 6, 3, 46,
          notes = "SRM column treated as higher management tier.")
add_cells("Lucocq_2024", "Reader 1 vs Reader 2, version 2019", "direct_reader_pair", 2, 7, 3, 49,
          notes = "SRM column treated as higher management tier.")
add_cells("Pitra_2018", "Reader 1 vs Reader 2, CT", "direct_reader_pair", 4, 3, 1, 21)
add_cells("Ragel_2016", "Observer 1 vs Observer 2, CT", "direct_reader_pair", 40, 6, 1, 4)
add_cells("Ragel_2016", "Observer 1 vs Observer 2, CEUS", "direct_reader_pair", 24, 18, 0, 9)
add_cells("Walmer_2023", "Reader 2 vs Reader 1", "direct_reader_pair", 24, 10, 3, 26,
          scale = "modified three-category CEUS proxy",
          notes = "Benign and indeterminate combined as lower tier; malignant as higher tier.")
add_cells("Walmer_2023", "Reader 3 vs Reader 1", "direct_reader_pair", 24, 9, 3, 27,
          scale = "modified three-category CEUS proxy")
add_cells("Walmer_2023", "Reader 4 vs Reader 1", "direct_reader_pair", 22, 9, 5, 27,
          scale = "modified three-category CEUS proxy")
add_cells("Walmer_2023", "Reader 3 vs Reader 2", "direct_reader_pair", 27, 6, 7, 23,
          scale = "modified three-category CEUS proxy")
add_cells("Walmer_2023", "Reader 4 vs Reader 2", "direct_reader_pair", 28, 3, 6, 26,
          scale = "modified three-category CEUS proxy")
add_cells("Walmer_2023", "Reader 4 vs Reader 3", "direct_reader_pair", 26, 5, 7, 25,
          scale = "modified three-category CEUS proxy")
add_cells("Yenice_2020", "Reader 1 vs Reader 2, CT", "direct_reader_pair", 58, 5, 5, 34)
add_cells("Yenice_2020", "Reader 1 vs Reader 2, MRI", "direct_reader_pair", 50, 7, 7, 38)

# Graumann: rows are gold-standard categories and columns are reader categories.
graumann <- list(
  "Reader A, session 1 vs gold standard" = matrix(c(
    18,0,0,1,0, 6,16,3,2,0, 1,3,7,1,1, 0,4,4,11,0, 0,0,0,5,17), 5, byrow = TRUE),
  "Reader A, session 2 vs gold standard" = matrix(c(
    18,0,0,1,0, 6,16,3,2,0, 1,3,6,2,1, 0,4,6,9,0, 0,0,1,4,17), 5, byrow = TRUE),
  "Reader B, session 1 vs gold standard" = matrix(c(
    18,1,0,0,0, 1,24,1,1,0, 0,0,12,1,0, 0,0,1,14,4, 0,0,0,0,22), 5, byrow = TRUE),
  "Reader B, session 2 vs gold standard" = matrix(c(
    18,1,0,0,0, 1,24,2,0,0, 0,0,9,4,0, 0,0,2,13,4, 0,0,0,0,22), 5, byrow = TRUE),
  "Reader C, session 1 vs gold standard" = matrix(c(
    19,0,0,0,0, 0,25,2,0,0, 0,3,10,0,0, 0,0,1,18,0, 0,0,0,0,22), 5, byrow = TRUE),
  "Reader C, session 2 vs gold standard" = matrix(c(
    19,0,0,0,0, 0,24,3,0,0, 0,3,10,0,0, 0,0,1,18,0, 0,0,0,0,22), 5, byrow = TRUE)
)
for (nm in names(graumann)) add_matrix("Graumann_2015", nm, "reader_reference", graumann[[nm]], 1:3)

# He: Tables 3-6 compare a reader with consensus; Tables 7-8 compare readers directly.
he <- list(
  "Low-experience reader 1 vs consensus" = list("reader_reference", matrix(c(
    122,22,0,0,0, 20,86,50,3,0, 0,13,111,21,0, 0,0,30,75,17, 0,0,0,6,27), 5, byrow = TRUE)),
  "Low-experience reader 2 vs consensus" = list("reader_reference", matrix(c(
    125,16,0,0,0, 17,65,63,3,0, 0,40,94,21,0, 0,0,34,71,4, 0,0,0,10,40), 5, byrow = TRUE)),
  "High-experience reader 3 vs consensus" = list("reader_reference", matrix(c(
    129,10,0,0,0, 13,101,9,0,0, 0,10,140,3,0, 0,0,42,98,0, 0,0,0,4,44), 5, byrow = TRUE)),
  "High-experience reader 4 vs consensus" = list("reader_reference", matrix(c(
    132,21,0,0,0, 10,85,12,0,0, 0,15,152,3,0, 0,0,27,100,4, 0,0,0,2,40), 5, byrow = TRUE)),
  "Low-experience reader 1 vs reader 2" = list("direct_reader_pair", matrix(c(
    108,25,11,0,0, 29,78,44,8,0, 4,39,72,29,1, 0,6,28,67,21, 0,0,0,5,28), 5, byrow = TRUE)),
  "High-experience reader 3 vs reader 4" = list("direct_reader_pair", matrix(c(
    127,11,1,0,0, 24,86,13,0,0, 2,10,132,9,0, 0,0,24,116,0, 0,0,0,6,42), 5, byrow = TRUE))
)
for (nm in names(he)) add_matrix("He_2026", nm, he[[nm]][[1]], he[[nm]][[2]], 1:3)

# Munch: rows are gold-standard categories and columns are reader categories.
munch <- list(
  "Expert reader 1 vs gold standard" = matrix(c(
    22,8,0,0,0, 0,67,4,0,1, 0,2,21,5,0, 0,0,0,14,2, 0,0,0,5,17), 5, byrow = TRUE),
  "Expert reader 2 vs gold standard" = matrix(c(
    27,2,1,0,0, 5,61,6,0,0, 0,1,25,2,0, 0,0,3,12,1, 0,0,0,1,21), 5, byrow = TRUE),
  "Intermediate reader 1 vs gold standard" = matrix(c(
    12,17,1,0,0, 9,46,11,4,2, 0,9,13,4,2, 0,1,5,10,0, 0,4,7,4,7), 5, byrow = TRUE),
  "Intermediate reader 2 vs gold standard" = matrix(c(
    17,11,2,0,0, 11,41,16,0,4, 1,10,12,3,2, 0,0,8,7,1, 0,4,8,6,4), 5, byrow = TRUE),
  "Beginner reader 2 vs gold standard" = matrix(c(
    18,9,0,3,0, 9,29,16,10,8, 0,4,15,6,3, 0,1,2,10,3, 0,2,0,5,15), 5, byrow = TRUE)
)
for (nm in names(munch)) add_matrix("Munch_2024", nm, "reader_reference", munch[[nm]], 1:3)
add_matrix("Munch_2024", "Beginner reader 1 vs beginner reader 2", "direct_reader_pair", matrix(c(
  3,1,0,0,0, 21,27,7,6,2, 3,13,22,15,11, 0,4,3,13,6, 0,0,1,0,10
), 5, byrow = TRUE), 1:3, scale = "modified EFSUMB CEUS Bosniak")
add_matrix("Munch_2024", "Intermediate reader 1 vs intermediate reader 2", "direct_reader_pair", matrix(c(
  15,5,1,0,0, 14,49,14,0,0, 0,12,23,2,0, 0,0,8,9,5, 0,0,0,5,6
), 5, byrow = TRUE), 1:3, scale = "modified EFSUMB CEUS Bosniak")
add_matrix("Munch_2024", "Expert reader 1 vs expert reader 2", "direct_reader_pair", matrix(c(
  20,2,0,0,0, 11,59,7,0,0, 0,2,22,1,0, 0,1,6,12,5, 1,0,0,2,17
), 5, byrow = TRUE), 1:3, scale = "modified EFSUMB CEUS Bosniak")
add_excluded(
  "Munch_2024", "Beginner reader 1 vs gold standard", "reader_reference", 168, "Bosniak",
  "Published cross-tabulation was internally inconsistent: the 25 printed cells and printed column totals summed to 177, whereas the printed grand total and study sample were 168.",
  "The BII row cells summed to 81 although its printed row total was 72; no cell was altered or imputed."
)

# Rosenkrantz: individual lesion assignments at 1.5 T and 3 T.
rosen <- data.frame(
  r1_15 = c("IIF","II","II","IV","II","II","II","II","II","IIF","II","III","IIF","IIF","IIF","II","IIF","II","IIF","III","II","II","III","I","IIF","IIF","II","II","IIF","II","II","III","IIF"),
  r1_30 = c("IV","IIF","II","IV","II","II","II","II","II","III","II","IIF","III","IIF","IIF","II","IIF","II","III","IV","II","II","IIF","II","II","IIF","IIF","IIF","III","II","II","III","II"),
  r2_15 = c("IIF","IIF","II","III","II","IIF","II","II","II","IIF","II","IIF","II","IIF","IIF","IIF","IIF","II","II","IV","II","IIF","IIF","II","II","III","IIF","IIF","IIF","II","II","IIF","II"),
  r2_30 = c("III","IIF","II","III","II","IIF","II","II","II","IIF","II","IIF","IIF","IIF","IIF","IIF","IIF","II","II","IV","IIF","IIF","IIF","II","II","III","IIF","IIF","IIF","II","II","IIF","IIF")
)
to_high <- function(x) x %in% c("III", "IV")
binary_cells <- function(a, b) {
  ah <- to_high(a); bh <- to_high(b)
  c(ll = sum(!ah & !bh), lh = sum(!ah & bh), hl = sum(ah & !bh), hh = sum(ah & bh))
}
cc <- binary_cells(rosen$r1_15, rosen$r2_15)
add_cells("Rosenkrantz_2014", "Reader 1 vs Reader 2, 1.5 T", "direct_reader_pair", cc[1], cc[2], cc[3], cc[4])
cc <- binary_cells(rosen$r1_30, rosen$r2_30)
add_cells("Rosenkrantz_2014", "Reader 1 vs Reader 2, 3 T", "direct_reader_pair", cc[1], cc[2], cc[3], cc[4])

add_excluded(
  "Pitra_2018", "Reader 1 vs Reader 2, MRI", "direct_reader_pair", 20, "Bosniak",
  "Binary kappa not estimable because all 20 valid paired classifications were in the higher management tier.",
  "Nine additional MRI table entries were non-Bosniak or missing."
)

cells <- do.call(rbind, records)
cells$effect_id <- sprintf("MTIR%03d", seq_len(nrow(cells)))
cells <- cells[, c("effect_id", setdiff(names(cells), "effect_id"))]

calc_stats <- function(ll, lh, hl, hh) {
  n <- ll + lh + hl + hh
  po <- (ll + hh) / n
  rlow <- ll + lh
  rhigh <- hl + hh
  clow <- ll + hl
  chigh <- lh + hh
  pe_k <- (rlow * clow + rhigh * chigh) / n^2
  kappa <- if (abs(1 - pe_k) < 1e-12) NA_real_ else (po - pe_k) / (1 - pe_k)
  p_low <- (rlow + clow) / (2 * n)
  pe_ac1 <- 2 * p_low * (1 - p_low)
  ac1 <- if (abs(1 - pe_ac1) < 1e-12) NA_real_ else (po - pe_ac1) / (1 - pe_ac1)
  c(observed_agreement = po, kappa = kappa, ac1 = ac1)
}

stats <- matrix(NA_real_, nrow(cells), 3L, dimnames = list(NULL, c("observed_agreement", "kappa", "ac1")))
for (i in which(cells$included == 1L)) {
  stats[i, ] <- calc_stats(cells$lower_lower[i], cells$lower_higher[i], cells$higher_lower[i], cells$higher_higher[i])
}
cells <- cbind(cells, as.data.frame(stats))
cells$yi_kappa <- ifelse(cells$included == 1L, atanh(pmax(pmin(cells$kappa, 0.995), -0.995)), NA_real_)
cells$yi_ac1 <- ifelse(cells$included == 1L, atanh(pmax(pmin(cells$ac1, 0.995), -0.995)), NA_real_)
cells$vi <- ifelse(cells$included == 1L & cells$n > 3, 1 / (cells$n - 3), NA_real_)

stopifnot(sum(cells$included) == 40L)
stopifnot(length(unique(cells$study_id[cells$included == 1L])) == 14L)
stopifnot(all(abs(rowSums(cells[cells$included == 1L, c("lower_lower", "lower_higher", "higher_lower", "higher_higher")]) -
                  cells$n[cells$included == 1L]) < 1e-10))

write.csv(cells, file.path(out_dir, "interreader_management_tier_cells.csv"), row.names = FALSE, na = "")

safe_tcrit <- function(df) if (is.finite(df) && df > 0) qt(0.975, df) else qnorm(0.975)

run_rve <- function(d, yi_col, label) {
  d <- d[d$included == 1L & is.finite(d[[yi_col]]) & is.finite(d$vi), ]
  model <- rma.mv(
    yi = d[[yi_col]],
    V = d$vi,
    random = list(~ 1 | study_id, ~ 1 | effect_id),
    data = d,
    method = "REML"
  )
  robust <- coef_test(model, vcov = "CR2", cluster = d$study_id)
  est_z <- as.numeric(model$b[1])
  se <- robust$SE[1]
  df <- robust$df_Satt[1]
  crit <- safe_tcrit(df)
  ci_z <- est_z + c(-1, 1) * crit * se
  pi_z <- est_z + c(-1, 1) * crit * sqrt(se^2 + sum(model$sigma2))
  data.frame(
    analysis = label,
    effects = nrow(d),
    studies = length(unique(d$study_id)),
    estimate = tanh(est_z),
    ci_lower = tanh(ci_z[1]),
    ci_upper = tanh(ci_z[2]),
    pi_lower = tanh(pi_z[1]),
    pi_upper = tanh(pi_z[2]),
    robust_se_z = se,
    satterthwaite_df = df,
    p_value = robust$p_Satt[1],
    tau2_study = model$sigma2[1],
    tau2_effect = model$sigma2[2]
  )
}

included <- cells[cells$included == 1L, ]
direct <- included[included$source_type == "direct_reader_pair", ]
bosniak_only <- included[included$scale == "Bosniak", ]
direct_bosniak <- included[included$source_type == "direct_reader_pair" & included$scale == "Bosniak", ]

summary_rows <- rbind(
  run_rve(direct_bosniak, "yi_kappa", "Direct reader pairs, Bosniak only: kappa"),
  run_rve(direct_bosniak, "yi_ac1", "Direct reader pairs, Bosniak only: AC1"),
  run_rve(included, "yi_kappa", "All reconstructed comparisons: kappa"),
  run_rve(included, "yi_ac1", "All reconstructed comparisons: AC1"),
  run_rve(direct, "yi_kappa", "Direct reader pairs only: kappa"),
  run_rve(direct, "yi_ac1", "Direct reader pairs only: AC1"),
  run_rve(bosniak_only, "yi_kappa", "Bosniak classifications only: kappa"),
  run_rve(bosniak_only, "yi_ac1", "Bosniak classifications only: AC1")
)
write.csv(summary_rows, file.path(out_dir, "interreader_management_tier_summary.csv"), row.names = FALSE)

descriptive <- do.call(rbind, lapply(list(
  "Direct reader pairs, Bosniak only" = direct_bosniak,
  "Direct reader pairs, all scales" = direct,
  "All reconstructed comparisons" = included
), function(d) {
  data.frame(
    analysis = names(d)[1],
    effects = nrow(d),
    studies = length(unique(d$study_id)),
    observed_agreement_min = min(d$observed_agreement),
    observed_agreement_median = median(d$observed_agreement),
    observed_agreement_max = max(d$observed_agreement),
    kappa_min = min(d$kappa),
    kappa_median = median(d$kappa),
    kappa_max = max(d$kappa),
    ac1_min = min(d$ac1),
    ac1_median = median(d$ac1),
    ac1_max = max(d$ac1)
  )
}))
# Preserve the list names explicitly; names(d) above refers to columns.
descriptive$analysis <- c(
  "Direct reader pairs, Bosniak only",
  "Direct reader pairs, all scales",
  "All reconstructed comparisons"
)
write.csv(descriptive, file.path(out_dir, "interreader_management_tier_descriptive.csv"), row.names = FALSE)

accounting <- data.frame(
  item = c(
    "has_raw_table inter-reader comparison rows",
    "duplicate metric or summary records",
    "unique reconstructable binary tables",
    "estimable binary tables",
    "included studies",
    "direct reader-pair tables",
    "reader-reference tables",
    "clinical-reassessment tables",
    "excluded binary tables"
  ),
  value = c(51L, 9L, nrow(cells), sum(cells$included), length(unique(included$study_id)),
            sum(included$source_type == "direct_reader_pair"),
            sum(included$source_type == "reader_reference"),
            sum(included$source_type == "clinical_reassessment"),
            sum(cells$included == 0L)),
  notes = c(
    "Rows in comparisons.csv; includes repeated metrics and summary rows.",
    "Kim 1, Lucocq 1, Ragel 1, and Walmer 6 records repeated an underlying raw table.",
    "One row per unique reconstructable 2x2 management-tier table, including excluded tables.",
    "Unique tables contributing a finite binary kappa and AC1.",
    "Studies contributing at least one estimable binary table.",
    "Direct reader-reader comparisons.",
    "Individual reader compared with consensus or gold standard.",
    "Initial report compared with specialist or multidisciplinary reassessment.",
    "Pitra 2018 MRI had degenerate binary margins; Munch 2024 beginner-reader-1 cells contradicted the printed totals."
  )
)
write.csv(accounting, file.path(out_dir, "interreader_management_tier_accounting.csv"), row.names = FALSE)

print(summary_rows)
print(accounting)
