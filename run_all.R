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
  "meta_analysis/01_main_analysis.R",
  "meta_analysis/02_additional_sensitivity.R",
  "meta_analysis/03_publication_bias.R",
  "meta_analysis/04_forest_plots.R",
  "meta_analysis/04b_forest_compact.R",
  "revision/analysis/revision_analyses.R"
)

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

message("\nAll analyses completed successfully.")
