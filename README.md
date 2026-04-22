# Predictors of High-Cost Healthcare Users in the US

**An Extension of Maynou et al. (2023) using MEPS 2022 Data**

This repository contains the complete analytical pipeline to predict high-cost healthcare users in the United States, replicating the methodological framework of Maynou, Street & García-Altés (2023) and extending it to examine the modulating role of health insurance.

## 📌 Project Overview

Identifying high-cost patients is critical for health system sustainability. While Maynou et al. (2023) demonstrated that clinical complexity (morbidity) dominates demographic factors in predicting high-cost status in Catalunya's universal system, this project tests whether that relationship holds in the fragmented US healthcare system.

**Key Extension:** We hypothesize that in the US, the morbidity-cost gradient is significantly modulated by insurance coverage type (Private, Medicare, Medicaid, Uninsured).

### Data Source
- **Dataset:** Medical Expenditure Panel Survey (MEPS) 2022 Full-Year Consolidated File (HC-243)
- **Sample:** Adults aged 18+
- **Outcome:** Total annual healthcare expenditure (`TOTEXP22`), operationalized as Top 5% and Top 1% thresholds.

## 📂 Repository Structure

```text
meps-high-cost-users/
├── R/                          # Analytical scripts (run in order)
│   ├── 01_download.R           # Downloads MEPS HC-243 and filters 18+
│   ├── 02_clean.R              # Variable construction and survey design
│   ├── 03_descriptive.R        # Summary statistics and zero-expenditure analysis
│   ├── 04_models_logistic.R    # Quasibinomial models for Top 5% / Top 1%
│   ├── 05_models_glm.R         # Gamma GLM (log link) for continuous expenditure
│   ├── 06_figures.R            # Data visualization (ggplot2)
│   └── 07_report.R             # Generates RMarkdown summary report
├── data/
│   ├── raw/                    # Raw .dta files (git-ignored)
│   └── processed/              # Cleaned .rds files and survey design objects
├── output/
│   ├── figures/                # Generated plots (.png, .pdf)
│   ├── models/                 # Saved model objects and AMEs
│   └── tables/                 # HTML tables and CSV results
├── docs/                       # Final analytical summary report
├── _run_all.R                  # Master script to execute the full pipeline
├── .gitignore                  # Standard R gitignore
└── README.md                   # This file
```

## 🚀 How to Run

1. **Prerequisites:** Ensure you have R installed along with the following packages:
   `haven`, `survey`, `dplyr`, `ggplot2`, `broom`, `marginaleffects`, `gtsummary`, `flextable`, `tidyr`, `scales`, `rmarkdown`.

2. **Execution:**
   Open the project in RStudio or your preferred IDE, set the working directory to the repository root, and run the master script:
   ```R
   source("_run_all.R")
   ```
   Alternatively, run the scripts in the `R/` folder sequentially from `01` to `07`.

## 📊 Methodology Highlights

1. **Survey Design:** All analyses incorporate MEPS complex survey weights (`PERWT22F`), strata (`VARSTR`), and PSUs (`VARPSU`) using the `survey` package.
2. **Thresholds:** Top 5% and Top 1% cutoffs are calculated using survey-weighted quantiles (`svyquantile`).
3. **Modeling:** 
   - Binary outcomes (Top 5/1%): Quasibinomial logistic regression.
   - Continuous outcomes: Gamma GLM with log link, strictly restricted to `total_exp > 0`.
4. **Interpretation:** Results are presented as Average Marginal Effects (AMEs) using the `marginaleffects` package, providing intuitive probability and dollar-scale interpretations.

## 📝 References
- Maynou, L., Street, A., & García-Altés, A. (2023). Profiling high-cost patients: A population-based study. *Health Policy*.
- Agency for Healthcare Research and Quality (AHRQ). Medical Expenditure Panel Survey (MEPS) 2022.
