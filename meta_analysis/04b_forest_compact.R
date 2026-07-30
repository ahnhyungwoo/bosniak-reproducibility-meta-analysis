## ============================================================
## Bosniak SRMA — Compact Forest Plots (Manuscript Main)
##
## Inter-reader: subgrouped by modality
## Inter-modality: subgrouped by modality pair
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
      tolower(sub("^.*_\\d{4}_", "", study_id)), ""),
    label = ifelse(suffix == "", paste0(author, " (", year, ")"),
      paste0(author, " (", year, suffix, ")")),
    mod_group = case_when(
      std_modality_1 == "CT" ~ "CT",
      std_modality_1 == "MRI" ~ "MRI",
      std_modality_1 %in% c("CEUS","US","SMI") ~ "CEUS/US",
      TRUE ~ "CT/MRI"
    ),
    mod_pair = case_when(
      (std_modality_1 %in% c("CT","DECT") & std_modality_2 == "MRI") |
        (std_modality_1 == "MRI" & std_modality_2 %in% c("CT","DECT")) ~ "CT vs MRI",
      (std_modality_1 %in% c("CT","CT_MRI") & std_modality_2 %in% c("CEUS","US")) |
        (std_modality_1 %in% c("CEUS","US") & std_modality_2 %in% c("CT","CT_MRI")) ~ "CT vs CEUS/US",
      (std_modality_1 == "MRI" & std_modality_2 %in% c("CEUS","US")) |
        (std_modality_1 %in% c("CEUS","US") & std_modality_2 == "MRI") ~ "MRI vs CEUS/US",
      TRUE ~ "Other"
    )
  ) %>% filter(!is.na(vi) & vi > 0)
}

opes_all <- comp %>%
  filter(is.na(exclude)|exclude=="", opes_include==1, !is.na(as.numeric(kappa))) %>%
  prep_opes()

# RVE estimates
dat_full <- comp %>%
  filter(is.na(exclude)|exclude=="",
         !is.na(as.numeric(kappa))) %>%
  mutate(
    stat = as.numeric(kappa),
    stat_clamped = pmin(pmax(stat,-0.995),0.995),
    yi = atanh(stat_clamped), nc = as.numeric(n_cysts),
    ci_lo=as.numeric(ci_lower), ci_hi=as.numeric(ci_upper), se_raw=as.numeric(kappa_se),
    vi_ci=ifelse(!is.na(ci_lo)&!is.na(ci_hi),((atanh(pmin(pmax(ci_hi,-0.995),0.995))-atanh(pmin(pmax(ci_lo,-0.995),0.995)))/(2*qnorm(0.975)))^2,NA_real_),
    vi_se=ifelse(!is.na(se_raw),(se_raw/(1-stat_clamped^2))^2,NA_real_),
    vi_n=ifelse(!is.na(nc)&nc>3,1/(nc-3),NA_real_),
    vi=coalesce(vi_ci,vi_se,vi_n)
  ) %>% filter(!is.na(stat)&!is.na(vi)&vi>0)

# Diamond estimates: true RVE+CR2 from all kappa effects
library(clubSandwich)
get_rve_pooled <- function(stratum) {
  sub <- dat_full %>% filter(analytic_stratum==stratum)
  if (nrow(sub) < 2) return(list(kappa=NA, ci_lo=NA, ci_hi=NA, k=nrow(sub), clusters=length(unique(sub$dep_id))))
  mod <- rma.mv(yi, vi, random=list(~1|dep_id, ~1|comparison_id), data=sub, method="REML")
  cr2 <- coef_test(mod, vcov="CR2", cluster=sub$dep_id)
  t_crit <- qt(0.975, cr2$df_Satt)
  ci_lo_z <- mod$b[1] - t_crit * cr2$SE
  ci_hi_z <- mod$b[1] + t_crit * cr2$SE
  list(kappa=tanh(mod$b[1]),
    ci_lo=tanh(ci_lo_z), ci_hi=tanh(ci_hi_z),
    k=nrow(sub), clusters=length(unique(sub$dep_id)))
}
rve_ir <- get_rve_pooled("inter-reader")
rve_im <- get_rve_pooled("inter-modality")


# ============================================================
# SUBGROUP FOREST PLOT FUNCTION
# ============================================================

make_subgroup_forest <- function(opes_data, group_var, group_order, group_labels,
                                 stratum_label, rve_est,
                                 filename, fig_width = 7.5,
                                 show_n = TRUE, show_group_col = FALSE,
                                 cex_val = 0.80, row_gap = 3,
                                 footnote = NULL) {

  # Sort by group then alphabetically
  opes_data$grp <- factor(opes_data[[group_var]], levels = group_order)
  opes_data <- opes_data %>% arrange(grp, label) %>% filter(!is.na(grp))
  k <- nrow(opes_data)

  # Fit overall RE model
  mod <- rma(yi, vi, data = opes_data, method = "REML", slab = opes_data$label)

  # Calculate row positions (bottom-up, with gaps between groups)
  groups <- group_order[group_order %in% opes_data$grp]
  n_groups <- length(groups)

  # Build row assignments (A at top = highest row within group)
  # row_gap = 3: 1 row for subtotal diamond + 2 rows visual gap
  row_pos <- numeric(k)
  group_header_rows <- numeric(n_groups)
  subtotal_rows <- numeric(n_groups)
  current_row <- 1
  for (g_idx in rev(seq_along(groups))) {
    g <- groups[g_idx]
    g_mask <- opes_data$grp == g
    n_in_group <- sum(g_mask)
    which_rows <- which(g_mask)
    # Reverse so first alphabetically gets highest row (top of group)
    row_pos[which_rows] <- (current_row + n_in_group - 1):current_row
    subtotal_rows[g_idx] <- current_row - 1
    group_header_rows[g_idx] <- current_row + n_in_group + 0.5
    current_row <- current_row + n_in_group + row_gap
  }

  # Pooled estimates at bottom
  row_opes <- -2.5
  row_rve <- -4

  # Total figure height
  max_row <- max(row_pos) + 3  # extra for top header
  fig_height <- max(max_row * 0.28 + 3.5, 6)

  # ilab: just n
  ilab_data <- NULL; ilab_pos <- NULL; ilab_head <- NULL
  if (show_n) {
    ilab_data <- cbind(ilab_data, opes_data$nc)
    ilab_pos <- c(ilab_pos, -0.60)
    ilab_head <- c(ilab_head, "n")
  }

  # Output both PDF and PNG
  for (ext in c("pdf", "png")) {
    fname <- sub("\\.pdf$", paste0(".", ext), filename)
    if (ext == "pdf") {
      pdf(fname, width = fig_width, height = fig_height)
    } else {
      png(fname, width = fig_width * 250, height = fig_height * 250, res = 250)
    }

    par(mar = c(4, 0, 1.5, 1), cex = cex_val, bg = "white")

    forest(mod, transf = tanh,
           header = c("Study", "\u03ba [95% CI]"),
           xlim = c(-1.4, 1.45),
           alim = c(-0.5, 1.0),
           at = seq(-0.4, 1.0, 0.2),
           refline = NA,
           rows = row_pos,
           ylim = c(row_rve - 1.5, max(row_pos) + 3),
           xlab = "",
           mlab = "",
           addfit = FALSE,
           ilab = ilab_data,
           ilab.xpos = ilab_pos,
           ilab.pos = 2,
           cex = cex_val,
           efac = 1.2,
           digits = 2)

    # Column header for n
    if (show_n) {
      text(ilab_pos[1], max(row_pos) + 2, "n", pos = 2, font = 2, cex = cex_val)
    }

    # Subgroup headers
    for (g_idx in seq_along(groups)) {
      g <- groups[g_idx]
      g_label <- group_labels[group_order == g]
      g_mask <- opes_data$grp == g
      n_g <- sum(g_mask)

      # Header text (bold, left-aligned)
      text(-1.4, group_header_rows[g_idx], pos = 4,
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

    # Overall OPES pooled diamond
    addpoly(mod$b[1], ci.lb = mod$ci.lb, ci.ub = mod$ci.ub,
            rows = row_opes, transf = tanh,
            mlab = sprintf("Overall OPES (k = %d)", k),
            col = "black", border = "black",
            cex = cex_val, efac = 1.3)

    # RVE+CR2 diamond
    addpoly(atanh(rve_est$kappa),
            ci.lb = atanh(rve_est$ci_lo), ci.ub = atanh(rve_est$ci_hi),
            rows = row_rve, transf = tanh,
            mlab = sprintf("Primary RVE+CR2 (k = %d, %d clusters)",
                           rve_est$k, rve_est$clusters),
            col = "firebrick4", border = "firebrick4",
            cex = cex_val, efac = 1.3)

    # Title
    title(stratum_label, line = 0.3, cex.main = 1.1, font.main = 1)

    # Footnote
    if (!is.null(footnote)) {
      mtext(footnote, side = 1, line = 2.8, adj = 0, cex = cex_val * 0.8, font = 3)
    }

    dev.off()
  }

  cat(sprintf("Saved: %s + .png (k=%d, %d groups, %.1f x %.1f in)\n",
              filename, k, n_groups, fig_width, fig_height))
}


# ============================================================
# FIGURE 1: Inter-reader, subgrouped by modality
# ============================================================

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("MANUSCRIPT FOREST PLOTS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

ir_opes <- opes_all %>% filter(analytic_stratum == "inter-reader")

# Check group sizes
cat("\nIR modality groups:\n")
print(table(ir_opes$mod_group))

make_subgroup_forest(
  ir_opes,
  group_var = "mod_group",
  group_order = c("CT", "MRI", "CEUS/US", "CT/MRI"),
  group_labels = c("CT", "MRI", "CEUS / US", "CT / MRI"),
  stratum_label = "Inter-reader agreement for Bosniak classification (\u03ba)",
  rve_est = rve_ir,
  filename = "fig1_inter_reader.pdf",
  fig_width = 7.5,
  cex_val = 0.80
)


# ============================================================
# FIGURE 2: Inter-modality, subgrouped by pair
# ============================================================

im_opes <- opes_all %>% filter(analytic_stratum == "inter-modality") %>%
  mutate(label = case_when(
    std_modality_1 == "CT_MRI" | std_modality_2 == "CT_MRI" ~ paste0(label, "*"),
    std_modality_1 == "DECT"   | std_modality_2 == "DECT"   ~ paste0(label, "\u2020"),
    TRUE ~ label
  ))

cat("\nIM modality pair groups:\n")
print(table(im_opes$mod_pair))

make_subgroup_forest(
  im_opes,
  group_var = "mod_pair",
  group_order = c("CT vs MRI", "CT vs CEUS/US", "MRI vs CEUS/US", "Other"),
  group_labels = c("CT vs MRI", "CT vs CEUS/US", "MRI vs CEUS/US", "Other"),
  stratum_label = "Inter-modality agreement for Bosniak classification (\u03ba)",
  rve_est = rve_im,
  filename = "fig2_inter_modality.pdf",
  fig_width = 7.5,
  cex_val = 0.80
)


cat("\nDone.\n")
