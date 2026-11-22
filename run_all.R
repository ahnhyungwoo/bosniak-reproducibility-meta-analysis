## Run all analyses from a clean repository checkout.

required_packages <- c("metafor", "clubSandwich", "dplyr", "readr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the required R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 0) stop("Run this file with Rscript.")
repo_dir <- dirname(
  normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
)

scripts <- c(
  "revision/analysis/revision_analyses.R",
  "revision/analysis/management_tier_interreader.R",
  "revision/analysis/management_tier_intermodality.R",
  "revision/analysis/variance_swap_sensitivity.R",
  "revision/analysis/reader_reference_sensitivity.R",
  "revision/analysis/R2_full_table_models.R",
  "revision/figures/r2/figure_data.R"
)
setwd(repo_dir)

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

for (relative_path in scripts) {
  script_path <- file.path(repo_dir, relative_path)
  message("\nRunning ", relative_path)
  status <- system2(rscript, shQuote(script_path))
  if (!identical(status, 0L)) {
    stop("Analysis failed in ", relative_path, " (exit status ", status, ").")
  }
}

stopifnot(file.copy("revision/figures/r2/data/publication_bias_summary.csv",
                   "revision/analysis/results/publication_bias_summary.csv", overwrite=TRUE))
message("\nR2 analyses complete. For final-size figures, run python revision/figures/r2/render_figures.py")
