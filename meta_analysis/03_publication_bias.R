## ============================================================
## Bosniak SRMA — Publication Bias / Small-Study Effects
##
## OPES-based: funnel plot, Egger's, Begg's, trim-and-fill
## ============================================================

library(metafor)
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

# --- OPES data prep ---
opes_all <- comp %>%
  filter(is.na(exclude) | exclude == "", opes_include == 1, !is.na(as.numeric(kappa))) %>%
  mutate(
    kappa_num = as.numeric(kappa),
    nc = as.numeric(n_cysts),
    ci_lo = as.numeric(ci_lower), ci_hi = as.numeric(ci_upper),
    se_raw = as.numeric(kappa_se),
    opes_val_str = ifelse(grepl("^(median_synthetic|reported_overall_restored):", opes_basis),
      sub("^[^:]+:", "", opes_basis), NA_character_),
    opes_val = coalesce(as.numeric(opes_val_str), kappa_num),
    opes_clamped = pmin(pmax(opes_val, -0.995), 0.995),
    yi = atanh(opes_clamped),
    vi_ci = ifelse(!is.na(ci_lo) & !is.na(ci_hi) & is.na(opes_val_str),
      ((atanh(pmin(pmax(ci_hi,-0.995),0.995)) - atanh(pmin(pmax(ci_lo,-0.995),0.995))) /
        (2*qnorm(0.975)))^2, NA_real_),
    vi_se = ifelse(!is.na(se_raw) & is.na(opes_val_str),
      (se_raw/(1-opes_clamped^2))^2, NA_real_),
    vi_n = ifelse(!is.na(nc) & nc > 3, 1/(nc-3), NA_real_),
    vi = coalesce(vi_ci, vi_se, vi_n),
    author = sub("_\\d{4}.*$", "", study_id),
    label = paste0(author, " (", year, ")")
  ) %>%
  filter(!is.na(vi) & vi > 0)


# ============================================================
# Helper: make one funnel plot (clean white background)
# ============================================================
make_funnel <- function(model, title, file_prefix, ylim_max = NULL, steps = 5) {
  # Determine clean y-axis limits
  if (is.null(ylim_max)) {
    se_max <- max(sqrt(model$vi), na.rm = TRUE)
    ylim_max <- ceiling(se_max * 10) / 10  # round up to nearest 0.1
  }
  for (ext in c("pdf", "png")) {
    fname <- sprintf("%s.%s", file_prefix, ext)
    if (ext == "pdf") {
      pdf(fname, width = 7, height = 6)
    } else {
      png(
        fname,
        width = 7 * DPI,
        height = 6 * DPI,
        res = DPI,
        type = "cairo",
        family = "Arial"
      )
    }
    par(mar = c(4.5, 4.5, 2, 1), bg = "white")
    funnel(model,
           xlab = expression(paste("Fisher ", italic(z), "-transformed ", kappa)),
           main = title,
           back = "white",
           shade = "white",
           ylim = c(ylim_max, 0),
           steps = steps)
    dev.off()
  }
  cat(sprintf("  Saved: %s.pdf/.png\n", file_prefix))
}


# ============================================================
# Helper: run full pub bias tests
# ============================================================
run_tests <- function(data, label) {
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat(sprintf("%s (k = %d)\n", label, nrow(data)))
  cat(paste(rep("=", 60), collapse = ""), "\n")

  mod <- rma(yi, vi, data = data, method = "REML", slab = data$label)

  cat("\n--- Egger's regression test ---\n")
  egger <- regtest(mod, model = "lm")
  cat(sprintf("  z = %.3f, p = %.4f\n", egger$zval, egger$pval))

  cat("\n--- Begg's rank correlation test ---\n")
  begg <- ranktest(mod)
  cat(sprintf("  tau = %.3f, p = %.4f\n", begg$tau, begg$pval))

  cat("\n--- Trim-and-fill ---\n")
  tf <- trimfill(mod)
  cat(sprintf("  Imputed studies: %d (side: %s)\n", tf$k0, tf$side))
  cat(sprintf("  Original:    kappa = %.3f (%.3f, %.3f)\n",
              tanh(mod$b[1]), tanh(mod$ci.lb), tanh(mod$ci.ub)))
  cat(sprintf("  After T&F:   kappa = %.3f (%.3f, %.3f)\n",
              tanh(tf$b[1]), tanh(tf$ci.lb), tanh(tf$ci.ub)))
  cat(sprintf("  Difference:  %.3f\n", abs(tanh(tf$b[1]) - tanh(mod$b[1]))))
  cat("\n")

  invisible(list(mod = mod, egger = egger, begg = begg, tf = tf))
}


# ============================================================
# Inter-reader OPES
# ============================================================
ir_opes <- opes_all %>% filter(analytic_stratum == "inter-reader")
ir_res <- run_tests(ir_opes, "Inter-reader agreement")

# IR: funnel only (T&F imputed 0 → no separate T&F plot needed)
make_funnel(ir_res$mod,
            "Inter-reader agreement",
            "fig_ir_funnel")


# ============================================================
# Inter-modality OPES
# ============================================================
im_opes <- opes_all %>% filter(analytic_stratum == "inter-modality")
im_res <- run_tests(im_opes, "Inter-modality concordance")

# IM: funnel + trim-and-fill funnel (1 study imputed)
make_funnel(im_res$mod,
            "Inter-modality concordance",
            "fig_im_funnel")

make_funnel(im_res$tf,
            "Inter-modality concordance (trim-and-fill)",
            "fig_im_funnel_tf")


# ============================================================
# SUMMARY
# ============================================================
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat(sprintf("\nInter-reader (k=%d):\n", nrow(ir_opes)))
cat(sprintf("  Egger: p=%.4f  |  Begg: p=%.4f  |  T&F imputed: %d  |  delta kappa: %.3f\n",
            ir_res$egger$pval, ir_res$begg$pval, ir_res$tf$k0,
            abs(tanh(ir_res$tf$b[1]) - tanh(ir_res$mod$b[1]))))

cat(sprintf("\nInter-modality (k=%d):\n", nrow(im_opes)))
cat(sprintf("  Egger: p=%.4f  |  Begg: p=%.4f  |  T&F imputed: %d  |  delta kappa: %.3f\n",
            im_res$egger$pval, im_res$begg$pval, im_res$tf$k0,
            abs(tanh(im_res$tf$b[1]) - tanh(im_res$mod$b[1]))))

cat("\nFigures generated:\n")
cat("  fig_ir_funnel.pdf/.png           — IR funnel (T&F imputed 0, text only)\n")
cat("  fig_im_funnel.pdf/.png           — IM funnel\n")
cat("  fig_im_funnel_tf.pdf/.png        — IM trim-and-fill funnel\n")

cat("\nDone.\n")
