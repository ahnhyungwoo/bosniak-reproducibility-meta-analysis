## ============================================================
## Bosniak SRMA — Supplement Forest Plots (Full Information)
##
## IR: subgrouped by modality, columns Wt/Ver/n
## IM: subgrouped by pair, columns Ver/n
## Intra: simple, columns Mod/Ver/n
## IV: simple, columns Mod/Ver/n
## ============================================================

library(metafor)
library(clubSandwich)
library(dplyr)
library(readr)

DPI <- 300

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 0) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
setwd(script_dir)

comp <- read_csv("comparisons.csv", show_col_types = FALSE)
stud <- read_csv("studies.csv", show_col_types = FALSE)
comp <- comp %>% left_join(stud %>% select(study_id, year), by = "study_id")

# ============================================================
# DATA PREP
# ============================================================

prep_opes <- function(data) {
  data %>% mutate(
    kappa_num = as.numeric(kappa),
    ci_lo = as.numeric(ci_lower), ci_hi = as.numeric(ci_upper),
    se_raw = as.numeric(kappa_se), nc = as.numeric(n_cysts),
    opes_val_str = ifelse(grepl("^(median_synthetic|reported_overall_restored):", opes_basis),
      sub("^[^:]+:", "", opes_basis), NA_character_),
    opes_val = coalesce(as.numeric(opes_val_str), kappa_num),
    opes_clamped = pmin(pmax(opes_val, -0.995), 0.995),
    yi = atanh(opes_clamped),
    vi_ci = ifelse(!is.na(ci_lo) & !is.na(ci_hi) & is.na(opes_val_str),
      ((atanh(pmin(pmax(ci_hi,-0.995),0.995))-atanh(pmin(pmax(ci_lo,-0.995),0.995)))/(2*qnorm(0.975)))^2, NA_real_),
    vi_se = ifelse(!is.na(se_raw) & is.na(opes_val_str), (se_raw/(1-opes_clamped^2))^2, NA_real_),
    vi_n = ifelse(!is.na(nc) & nc > 3, 1/(nc-3), NA_real_),
    vi = coalesce(vi_ci, vi_se, vi_n),
    author = sub("_\\d{4}.*$", "", study_id),
    suffix = ifelse(grepl("_\\d{4}_.+$", study_id),
      toupper(sub("^.*_\\d{4}_", "", study_id)), ""),
    label = ifelse(suffix == "", paste0(author, " (", year, ")"),
      paste0(author, " (", year, suffix, ")")),
    wt_abbr = case_when(
      weight_scheme == "unweighted" ~ "UW", weight_scheme == "linear" ~ "LW",
      weight_scheme == "quadratic" ~ "QW", weight_scheme == "Cicchetti" ~ "Cic",
      weight_scheme == "weighted_unspecified" ~ "W-NS",
      weight_scheme == "unspecified" ~ "NS", weight_scheme == "custom" ~ "Cust",
      TRUE ~ "\u2014"),
    mod_abbr = case_when(
      std_modality_1 == "CT" ~ "CT", std_modality_1 == "MRI" ~ "MRI",
      std_modality_1 == "CEUS" ~ "CEUS", std_modality_1 == "US" ~ "US",
      std_modality_1 == "CT_MRI" ~ "CT/MRI", std_modality_1 == "SMI" ~ "SMI",
      std_modality_1 == "DECT" ~ "DECT", std_modality_1 == "MV_flow" ~ "MVF",
      TRUE ~ "\u2014"),
    ver_abbr = case_when(
      std_version_1 == "original" ~ "Original", std_version_1 == "v2019" ~ "v2019",
      std_version_1 == "modified" ~ "Modified", std_version_1 == "EFSUMB" ~ "EFSUMB",
      std_version_1 == "binary" ~ "Binary",
      std_version_1 == "both" ~ "Original/v2019", std_version_1 == "unspecified" ~ "NS",
      TRUE ~ "\u2014"),
    mod_group = case_when(
      std_modality_1 == "CT" ~ "CT", std_modality_1 == "MRI" ~ "MRI",
      std_modality_1 %in% c("CEUS","US","SMI") ~ "CEUS/US",
      TRUE ~ "CT/MRI"),
    mod_pair = case_when(
      (std_modality_1 %in% c("CT","DECT") & std_modality_2 == "MRI") |
        (std_modality_1 == "MRI" & std_modality_2 %in% c("CT","DECT")) ~ "CT vs MRI",
      (std_modality_1 %in% c("CT","CT_MRI") & std_modality_2 %in% c("CEUS","US")) |
        (std_modality_1 %in% c("CEUS","US") & std_modality_2 %in% c("CT","CT_MRI")) ~ "CT vs CEUS/US",
      (std_modality_1 == "MRI" & std_modality_2 %in% c("CEUS","US")) |
        (std_modality_1 %in% c("CEUS","US") & std_modality_2 == "MRI") ~ "MRI vs CEUS/US",
      TRUE ~ "Other")
  ) %>% filter(!is.na(vi) & vi > 0)
}

opes_all <- comp %>%
  filter(is.na(exclude)|exclude=="", opes_include==1, !is.na(as.numeric(kappa))) %>%
  prep_opes()

# RVE estimates
dat_full <- comp %>%
  filter(is.na(exclude)|exclude=="") %>%
  mutate(
    stat = as.numeric(kappa), stat_clamped = pmin(pmax(stat,-0.995),0.995),
    yi = atanh(stat_clamped), nc = as.numeric(n_cysts),
    ci_lo=as.numeric(ci_lower), ci_hi=as.numeric(ci_upper), se_raw=as.numeric(kappa_se),
    vi_ci=ifelse(!is.na(ci_lo)&!is.na(ci_hi),((atanh(pmin(pmax(ci_hi,-0.995),0.995))-atanh(pmin(pmax(ci_lo,-0.995),0.995)))/(2*qnorm(0.975)))^2,NA_real_),
    vi_se=ifelse(!is.na(se_raw),(se_raw/(1-stat_clamped^2))^2,NA_real_),
    vi_n=ifelse(!is.na(nc)&nc>3,1/(nc-3),NA_real_),
    vi=coalesce(vi_ci,vi_se,vi_n)
  ) %>% filter(!is.na(stat)&!is.na(vi)&vi>0)

# Diamond estimates: use OPES (matches figure legend "one-per-study estimates")
get_opes_pooled <- function(stratum) {
  sub <- opes_all %>% filter(analytic_stratum==stratum)
  if (nrow(sub) < 2) return(list(kappa=NA, ci_lo=NA, ci_hi=NA, k=nrow(sub), clusters=nrow(sub)))
  mod <- rma(yi, vi, data=sub, method="REML")
  list(kappa=tanh(mod$b[1]),
    ci_lo=tanh(mod$ci.lb), ci_hi=tanh(mod$ci.ub),
    k=nrow(sub), clusters=nrow(sub))
}
rve_ir <- get_opes_pooled("inter-reader")
rve_im <- get_opes_pooled("inter-modality")


# ============================================================
# SUBGROUP FOREST (with info columns)
# ============================================================

make_subgroup_forest_supp <- function(opes_data, group_var, group_order, group_labels,
                                      stratum_label, rve_est,
                                      show_cols, col_positions,
                                      filename, fig_width = 9,
                                      xlim_range = c(-1.6, 1.5),
                                      cex_val = 0.80, row_gap = 3,
                                      footnote = NULL) {

  col_map <- list(
    wt  = list(header = "Wt",  field = "wt_abbr"),
    mod = list(header = "Modality", field = "mod_abbr"),
    ver = list(header = "Version", field = "ver_abbr"),
    n   = list(header = "n",   field = "nc")
  )

  opes_data$grp <- factor(opes_data[[group_var]], levels = group_order)
  opes_data <- opes_data %>% arrange(grp, label) %>% filter(!is.na(grp))
  k <- nrow(opes_data)

  mod <- rma(yi, vi, data = opes_data, method = "REML", slab = opes_data$label)

  # Build ilab matrix from sorted data
  ilab_mat <- NULL; ilab_headers <- character(0)
  for (col in show_cols) {
    ilab_mat <- cbind(ilab_mat, opes_data[[col_map[[col]]$field]])
    ilab_headers <- c(ilab_headers, col_map[[col]]$header)
  }

  # Row positions (A at top within each group)
  # row_gap = 3: 1 row for subtotal diamond + 2 rows visual gap
  groups <- group_order[group_order %in% opes_data$grp]
  n_groups <- length(groups)
  row_pos <- numeric(k)
  group_header_rows <- numeric(n_groups)
  subtotal_rows <- numeric(n_groups)
  current_row <- 1
  for (g_idx in rev(seq_along(groups))) {
    g <- groups[g_idx]
    g_mask <- opes_data$grp == g
    n_in_group <- sum(g_mask)
    which_rows <- which(g_mask)
    row_pos[which_rows] <- (current_row + n_in_group - 1):current_row
    subtotal_rows[g_idx] <- current_row - 1
    group_header_rows[g_idx] <- current_row + n_in_group + 0.5
    current_row <- current_row + n_in_group + row_gap
  }

  row_opes <- -2.5
  row_rve <- -4
  max_row <- max(row_pos) + 3
  fig_height <- max(max_row * 0.28 + 3.5, 6)

  for (ext in c("pdf", "png")) {
    fname <- sub("\\.pdf$", paste0(".", ext), filename)
    if (ext == "pdf") pdf(fname, width = fig_width, height = fig_height)
    else png(
      fname,
      width = fig_width * DPI,
      height = fig_height * DPI,
      res = DPI,
      type = "cairo",
      family = "Arial"
    )

    par(mar = c(4, 0, 1.5, 1), cex = cex_val, bg = "white")

    forest(mod, transf = tanh,
           header = c("Study", "\u03ba [95% CI]"),
           xlim = xlim_range,
           alim = c(-0.5, 1.0),
           at = seq(-0.4, 1.0, 0.2),
           refline = NA,
           rows = row_pos,
           ylim = c(row_rve - 1.5, max(row_pos) + 3),
           xlab = "",
           mlab = "",
           addfit = FALSE,
           ilab = ilab_mat,
           ilab.xpos = col_positions,
           ilab.pos = 2,
           cex = cex_val, efac = 1.2, digits = 2)

    # Column headers
    for (i in seq_along(ilab_headers))
      text(col_positions[i], max(row_pos) + 2, ilab_headers[i],
           pos = 2, font = 2, cex = cex_val)

    # Subgroup headers (bold)
    for (g_idx in seq_along(groups)) {
      g <- groups[g_idx]
      g_label <- group_labels[group_order == g]
      n_g <- sum(opes_data$grp == g)
      text(xlim_range[1], group_header_rows[g_idx], pos = 4,
           bquote(bold(.(paste0(g_label, " (k = ", n_g, ")")))),
           cex = cex_val * 1.05)
    }

    # Subgroup subtotal diamonds
    for (g_idx in seq_along(groups)) {
      g <- groups[g_idx]
      g_mask <- opes_data$grp == g
      g_data <- opes_data[g_mask, ]
      if (nrow(g_data) >= 2) {
        g_mod <- rma(yi, vi, data = g_data, method = "REML")
        addpoly(g_mod$b[1], ci.lb = g_mod$ci.lb, ci.ub = g_mod$ci.ub,
                rows = subtotal_rows[g_idx], transf = tanh,
                mlab = "Subtotal",
                col = "gray30", border = "gray30",
                cex = cex_val, efac = 0.8)
      }
    }

    # Separator line between subtotals and overall diamonds
    abline(h = -0.5, lty = 1, col = "gray50")

    # Overall OPES diamond
    addpoly(mod$b[1], ci.lb = mod$ci.lb, ci.ub = mod$ci.ub,
            rows = row_opes, transf = tanh,
            mlab = sprintf("Overall OPES (k = %d)", k),
            col = "black", border = "black", cex = cex_val, efac = 1.3)

    # RVE diamond
    if (!is.null(rve_est)) {
      addpoly(atanh(rve_est$kappa),
              ci.lb = atanh(rve_est$ci_lo), ci.ub = atanh(rve_est$ci_hi),
              rows = row_rve, transf = tanh,
              mlab = sprintf("Primary RVE+CR2 (k = %d, %d clusters)",
                             rve_est$k, rve_est$clusters),
              col = "firebrick4", border = "firebrick4", cex = cex_val, efac = 1.3)
    }

    title(stratum_label, line = 0.3, cex.main = 1.1, font.main = 1)

    if (!is.null(footnote)) {
      mtext(footnote, side = 1, line = 2.8, adj = 0, cex = cex_val * 0.8, font = 3)
    }

    dev.off()
  }
  cat(sprintf("Saved: %s + .png (k=%d, %d groups, %.1f x %.1f in)\n",
              filename, k, n_groups, fig_width, fig_height))
}


# ============================================================
# SIMPLE FOREST (no subgroups, for small strata)
# ============================================================

make_simple_forest_supp <- function(opes_data, stratum_label,
                                    show_cols, col_positions,
                                    filename, fig_width = 8,
                                    xlim_range = c(-1.5, 1.6),
                                    alim_range = c(-0.5, 1.0),
                                    at_seq = seq(-0.4, 1.0, 0.2),
                                    cex_val = 0.85) {

  col_map <- list(
    wt  = list(header = "Wt",  field = "wt_abbr"),
    mod = list(header = "Modality", field = "mod_abbr"),
    ver = list(header = "Version", field = "ver_abbr"),
    n   = list(header = "n",   field = "nc")
  )

  opes_data <- opes_data %>% arrange(label)
  k <- nrow(opes_data)
  if (k == 0) { cat("No data for", stratum_label, "\n"); return() }

  mod <- rma(yi, vi, data = opes_data, method = "REML", slab = opes_data$label)

  ilab_mat <- NULL; ilab_headers <- character(0)
  for (col in show_cols) {
    ilab_mat <- cbind(ilab_mat, opes_data[[col_map[[col]]$field]])
    ilab_headers <- c(ilab_headers, col_map[[col]]$header)
  }

  fig_height <- max(k * 0.45 + 3, 5)

  for (ext in c("pdf", "png")) {
    fname <- sub("\\.pdf$", paste0(".", ext), filename)
    if (ext == "pdf") pdf(fname, width = fig_width, height = fig_height)
    else png(
      fname,
      width = fig_width * DPI,
      height = fig_height * DPI,
      res = DPI,
      type = "cairo",
      family = "Arial"
    )

    par(mar = c(4, 0, 1.5, 1), cex = cex_val, bg = "white")

    forest(mod, transf = tanh,
           header = c("Study", "\u03ba [95% CI]"),
           xlim = xlim_range, alim = alim_range, at = at_seq,
           refline = NA, rows = k:1, ylim = c(-2, k + 3),
           xlab = "",
           mlab = sprintf("OPES random-effects (k = %d)", k),
           ilab = ilab_mat, ilab.xpos = col_positions, ilab.pos = 2,
           cex = cex_val, efac = 1.3, digits = 2)

    for (i in seq_along(ilab_headers))
      text(col_positions[i], k + 2, ilab_headers[i],
           pos = 2, font = 2, cex = cex_val)

    title(stratum_label, line = 0.3, cex.main = 1.1, font.main = 1)
    dev.off()
  }
  cat(sprintf("Saved: %s + .png (k=%d, %.1f x %.1f in)\n",
              filename, k, fig_width, fig_height))
}


# ============================================================
# GENERATE SUPPLEMENT PLOTS
# ============================================================

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("SUPPLEMENT FOREST PLOTS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# --- Inter-reader: subgrouped by modality, Wt/Ver/n ---
ir_opes <- opes_all %>% filter(analytic_stratum == "inter-reader")

make_subgroup_forest_supp(
  ir_opes,
  group_var = "mod_group",
  group_order = c("CT", "MRI", "CEUS/US", "CT/MRI"),
  group_labels = c("CT", "MRI", "CEUS / US", "CT / MRI"),
  stratum_label = "Inter-reader agreement for Bosniak classification (\u03ba)",
  rve_est = rve_ir,
  show_cols = c("wt", "ver", "n"),
  col_positions = c(-0.85, -0.70, -0.58),
  filename = "forest_inter_reader.pdf",
  fig_width = 9.5,
  xlim_range = c(-1.8, 1.5)
)

# --- Inter-modality: subgrouped by pair, Ver/n ---
im_opes <- opes_all %>% filter(analytic_stratum == "inter-modality") %>%
  mutate(label = case_when(
    std_modality_1 == "CT_MRI" | std_modality_2 == "CT_MRI" ~ paste0(label, "*"),
    std_modality_1 == "DECT"   | std_modality_2 == "DECT"   ~ paste0(label, "\u2020"),
    TRUE ~ label
  ))

make_subgroup_forest_supp(
  im_opes,
  group_var = "mod_pair",
  group_order = c("CT vs MRI", "CT vs CEUS/US", "MRI vs CEUS/US", "Other"),
  group_labels = c("CT vs MRI", "CT vs CEUS/US", "MRI vs CEUS/US", "Other"),
  stratum_label = "Inter-modality agreement for Bosniak classification (\u03ba)",
  rve_est = rve_im,
  show_cols = c("ver", "n"),
  col_positions = c(-0.70, -0.58),
  filename = "forest_inter_modality.pdf",
  fig_width = 9,
  xlim_range = c(-1.6, 1.5)
)

# --- Intra-reader: simple, Mod/n (Ver omitted for consistency with IV plot) ---
intra_opes <- opes_all %>% filter(analytic_stratum == "intra-reader")

make_simple_forest_supp(
  intra_opes,
  stratum_label = "Intra-reader agreement for Bosniak classification (\u03ba)",
  show_cols = c("mod", "n"),
  col_positions = c(-0.75, -0.58),
  filename = "forest_intra_reader.pdf",
  fig_width = 8.5,
  xlim_range = c(-1.6, 1.6),
  alim_range = c(0.0, 1.0),
  at_seq = seq(0.0, 1.0, 0.2)
)

# --- Inter-version: simple, Mod/n (Ver omitted: all are original vs v2019) ---
iv_opes <- opes_all %>% filter(analytic_stratum == "inter-version")

make_simple_forest_supp(
  iv_opes,
  stratum_label = "Inter-version agreement for Bosniak classification (\u03ba)",
  show_cols = c("mod", "n"),
  col_positions = c(-0.75, -0.58),
  filename = "forest_inter_version.pdf",
  fig_width = 8.5,
  xlim_range = c(-1.6, 1.6)
)

cat("\nDone.\n")
