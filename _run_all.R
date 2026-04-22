# _run_all.R
# MEPS High-Cost Patient Prediction Project
# Master script: runs the full pipeline in order.
# Run from the project root directory.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("R/01_download.R")
source("R/02_clean.R")
source("R/03_descriptive.R")
source("R/04_models_logistic.R")
source("R/05_models_glm.R")
source("R/06_figures.R")
source("R/07_report.R")

cat("\nPipeline complete. See output/ for results.\n")
