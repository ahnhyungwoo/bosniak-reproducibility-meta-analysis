## 동일 효과 집합에서 분산만 바꿔 1/(n-3) 근사의 영향을 직접 본다.
##
## 기존 'directly reported variance only' 민감도는 근사 분산을 쓴 효과를 통째로
## 제외한 부분집합 분석이므로, 차이에 연구 구성 변화가 섞인다. 여기서는 효과
## 집합을 고정하고 분산만 보고값에서 1/(n-3)로 바꿔 두 번 적합한다.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
project_dir <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/")

src <- readLines(file.path(project_dir, "revision", "analysis", "revision_analyses.R"))
stop_at <- grep("summary_parts <- list(", src, fixed = TRUE)[1]
eval(parse(text = paste(src[seq_len(stop_at - 1)], collapse = "\n")))

compare <- function(d, label) {
  # 두 모형이 같은 효과를 쓰도록 보고 분산과 n 근사가 모두 가능한 것만 남긴다.
  reported <- d %>% filter(variance_source %in% c("reported_CI", "reported_SE"))
  base <- reported %>% filter(!is.na(vi_z), vi_z > 0, !is.na(vi_z_n), vi_z_n > 0)
  dropped <- nrow(reported) - nrow(base)

  a <- run_rve(base, paste(label, "reported variance"),
               analysis_family = "variance_swap")
  b <- run_rve(base %>% mutate(vi_z = vi_z_n),
               paste(label, "1/(n-3) substituted"),
               analysis_family = "variance_swap")

  cat(sprintf("\n== %s ==\n", label))
  cat(sprintf("  고정 효과 집합 %d개 (n 결측으로 제외 %d개)\n", nrow(base), dropped))
  for (r in list(a, b)) {
    cat(sprintf("  %-24s k=%3d  kappa=%.4f  CI %.4f-%.4f  PI %.4f-%.4f  df=%.2f\n",
                sub(paste0(label, " "), "", r$analysis), r$effects,
                r$estimate, r$ci_lower, r$ci_upper, r$pi_lower, r$pi_upper,
                r$satterthwaite_df))
  }
  cat(sprintf("  차이  kappa %+.4f\n", b$estimate - a$estimate))
  bind_rows(a, b)
}

res <- bind_rows(compare(ir, "Inter-reader"), compare(im, "Inter-modality"))
write_csv(res, file.path(project_dir, "revision", "analysis", "results",
                         "variance_swap_sensitivity.csv"))
cat("\nresults/variance_swap_sensitivity.csv 저장\n")
