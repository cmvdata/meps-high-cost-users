# =============================================================================
# 06_twopart_profiles.R
#
# Purpose:
#   (1) Recompute weighted AMEs for GLM models G-B and G-C (closes Step 05).
#   (2) Estimate a participation model Pr(Y > 0 | X) with specification B.
#   (3) Build a 6-profile two-part expected-spending table combining:
#         - Pr(Y > 0 | X)   from the participation logit
#         - Pr(Top 5% | X)  from mod_B1
#         - Pr(Top 1% | X)  from mod_B2
#         - E[Y | Y > 0, X] from mod_GB (Gamma GLM, log link)
#         - E[Y | X] = Pr(Y > 0 | X) * E[Y | Y > 0, X]
#
# Inputs:
#   data/processed/meps_design.rds
#   output/models/logistic_models.rds   (list with B1, B2, C1, C2)
#   output/models/glm_models.rds        (list with GB, GC)
#
# Outputs:
#   output/tables/ame_GB_weighted.csv
#   output/tables/ame_GC_weighted.csv
#   output/tables/ame_participation.csv
#   output/tables/table_twopart_profiles.csv
#   output/tables/table_twopart_profiles.html
#   output/models/participation_model.rds
#
# Notes:
#   - Uncertainty for E[Y] is NOT propagated (point estimate only).
#     A methodological footnote is included in the HTML table.
#   - Participation model uses Specification B (no interaction), matching
#     mod_B1 and mod_GB for coherence.
# =============================================================================

suppressPackageStartupMessages({
  library(survey)
  library(marginaleffects)
  library(dplyr)
  library(tidyr)
  library(gt)
})

# -----------------------------------------------------------------------------
# 0. Load persisted objects
# -----------------------------------------------------------------------------
cat("Loading design and fitted models...\n")

design      <- readRDS("data/processed/meps_design.rds")
logit_list  <- readRDS("output/models/logistic_models.rds")
glm_list    <- readRDS("output/models/glm_models.rds")

mod_B1 <- logit_list$B1
mod_B2 <- logit_list$B2
mod_GB <- glm_list$GB
mod_GC <- glm_list$GC

# Recreate the positive-expenditure subdesign inside this script
design_pos <- subset(design, total_exp > 0)

# -----------------------------------------------------------------------------
# 1. Weighted AMEs for G-B and G-C (closes Step 05)
# -----------------------------------------------------------------------------
cat("\nRecomputing weighted AMEs for GLM models...\n")

ame_GB <- avg_slopes(
  mod_GB,
  wts = model.frame(mod_GB)$`(weights)`
)

ame_GC <- avg_slopes(
  mod_GC,
  variables = "n_chronic_conditions_cat",
  by        = "insurance_status",
  wts       = model.frame(mod_GC)$`(weights)`
)

write.csv(as.data.frame(ame_GB),
          "output/tables/ame_GB_weighted.csv", row.names = FALSE)
write.csv(as.data.frame(ame_GC),
          "output/tables/ame_GC_weighted.csv", row.names = FALSE)

cat("Weighted AMEs saved (G-B, G-C).\n")

# -----------------------------------------------------------------------------
# 2. Participation model: Pr(Y > 0 | X), specification B
# -----------------------------------------------------------------------------
cat("\nFitting participation model (Pr(Y > 0) | X)...\n")

# Build the binary participation indicator on the full design.
# Using update() so the variable lives inside the svydesign object.
design <- update(design, has_exp = as.numeric(total_exp > 0))

mod_P <- svyglm(
  has_exp ~ age_cat + sex + race_ethnicity + education +
            income_pct_fpl + insurance_status +
            n_chronic_conditions_cat + mental_health_cat + cancer_ever,
  design = design,
  family = quasibinomial()
)

# Convergence check (defensive: svyglm can return NULL or non-logical)
if (!isTRUE(mod_P$converged)) {
  warning("Participation model may not have converged. Inspect before using.")
} else {
  cat(sprintf("Participation model converged. n = %d\n",
              length(mod_P$y)))
}

# Deviance-based pseudo-R^2 for the participation model
null_dev_P  <- mod_P$null.deviance
res_dev_P   <- mod_P$deviance
pseudo_R2_P <- 1 - res_dev_P / null_dev_P
cat(sprintf("Participation pseudo-R^2 (deviance): %.4f\n", pseudo_R2_P))

# Weighted AMEs for the participation model
ame_P <- avg_slopes(
  mod_P,
  wts = model.frame(mod_P)$`(weights)`
)
write.csv(as.data.frame(ame_P),
          "output/tables/ame_participation.csv", row.names = FALSE)

# Persist the fitted model
saveRDS(mod_P, "output/models/participation_model.rds")
cat("Participation model saved to output/models/participation_model.rds\n")

# -----------------------------------------------------------------------------
# 3. Build 6-profile newdata with exact factor levels
# -----------------------------------------------------------------------------
cat("\nConstructing 6-profile newdata...\n")

# Helper: extract factor levels from the model frame to guarantee matching
get_levels <- function(mod, var) {
  mf <- model.frame(mod)
  if (is.factor(mf[[var]])) levels(mf[[var]]) else unique(mf[[var]])
}

# Pull levels from mod_GB (most complete spec with all B covariates)
lvl_age  <- get_levels(mod_GB, "age_cat")
lvl_sex  <- get_levels(mod_GB, "sex")
lvl_race <- get_levels(mod_GB, "race_ethnicity")
lvl_edu  <- get_levels(mod_GB, "education")
lvl_inc  <- get_levels(mod_GB, "income_pct_fpl")
lvl_ins  <- get_levels(mod_GB, "insurance_status")
lvl_ncc  <- get_levels(mod_GB, "n_chronic_conditions_cat")
lvl_mh   <- get_levels(mod_GB, "mental_health_cat")

# Defensive check: confirm the literal labels we plan to use exist
need_in <- function(x, pool, name) {
  if (!all(x %in% pool)) {
    stop(sprintf(
      "Level mismatch for %s. Requested: %s. Available: %s",
      name,
      paste(x, collapse = ", "),
      paste(pool, collapse = ", ")
    ))
  }
}
need_in("55-64",                  lvl_age,  "age_cat")
need_in("Female",                 lvl_sex,  "sex")
need_in("White Non-Hispanic",     lvl_race, "race_ethnicity")
need_in("College+",               lvl_edu,  "education")
need_in("Middle Income",          lvl_inc,  "income_pct_fpl")
need_in(c("Private",
          "Medicaid/Public Only",
          "Uninsured"),           lvl_ins,  "insurance_status")
need_in(c("0", "4+"),             lvl_ncc,  "n_chronic_conditions_cat")
need_in("No SPD",                 lvl_mh,   "mental_health_cat")

# Detect class of cancer_ever in the model frame; use matching type in newdata
cancer_class <- class(model.frame(mod_GB)$cancer_ever)[1]
cancer_no <- if (cancer_class == "factor") {
  lv <- levels(model.frame(mod_GB)$cancer_ever)
  if ("0"  %in% lv) "0" else if ("No" %in% lv) "No" else lv[1]
} else {
  0  # numeric
}
cat(sprintf("cancer_ever class in model frame: %s -> using %s\n",
            cancer_class, as.character(cancer_no)))

# Construct profile grid: 2 morbidity x 3 insurance = 6 rows
profiles <- expand.grid(
  n_chronic_conditions_cat = c("0", "4+"),
  insurance_status         = c("Private", "Medicaid/Public Only", "Uninsured"),
  age_cat                  = "55-64",
  sex                      = "Female",
  race_ethnicity           = "White Non-Hispanic",
  education                = "College+",
  income_pct_fpl           = "Middle Income",
  mental_health_cat        = "No SPD",
  stringsAsFactors         = FALSE
)
profiles$cancer_ever <- cancer_no

# Coerce to factors with the exact levels from the model frame
profiles$n_chronic_conditions_cat <- factor(profiles$n_chronic_conditions_cat,
                                            levels = lvl_ncc)
profiles$insurance_status         <- factor(profiles$insurance_status,
                                            levels = lvl_ins)
profiles$age_cat                  <- factor(profiles$age_cat,  levels = lvl_age)
profiles$sex                      <- factor(profiles$sex,      levels = lvl_sex)
profiles$race_ethnicity           <- factor(profiles$race_ethnicity,
                                            levels = lvl_race)
profiles$education                <- factor(profiles$education, levels = lvl_edu)
profiles$income_pct_fpl           <- factor(profiles$income_pct_fpl,
                                            levels = lvl_inc)
profiles$mental_health_cat        <- factor(profiles$mental_health_cat,
                                            levels = lvl_mh)
if (cancer_class == "factor") {
  profiles$cancer_ever <- factor(profiles$cancer_ever,
                                 levels = levels(model.frame(mod_GB)$cancer_ever))
}

cat("Profiles constructed:\n")
print(profiles[, c("insurance_status", "n_chronic_conditions_cat")])

# -----------------------------------------------------------------------------
# 4. Predictions from each model
# -----------------------------------------------------------------------------
cat("\nGenerating predictions for each model...\n")

pred_P  <- predictions(mod_P,  newdata = profiles, type = "response")
pred_B1 <- predictions(mod_B1, newdata = profiles, type = "response")
pred_B2 <- predictions(mod_B2, newdata = profiles, type = "response")
pred_GB <- predictions(mod_GB, newdata = profiles, type = "response")

# -----------------------------------------------------------------------------
# 5. Assemble the two-part table
# -----------------------------------------------------------------------------
cat("Assembling two-part expected-spending table...\n")

out <- profiles %>%
  mutate(
    pr_pos     = pred_P$estimate,
    pr_pos_se  = pred_P$std.error,
    pr_top5    = pred_B1$estimate,
    pr_top5_se = pred_B1$std.error,
    pr_top1    = pred_B2$estimate,
    pr_top1_se = pred_B2$std.error,
    e_y_pos    = pred_GB$estimate,
    e_y_pos_se = pred_GB$std.error,
    e_y        = pr_pos * e_y_pos    # two-part product; no SE by design
  ) %>%
  select(insurance_status, n_chronic_conditions_cat,
         pr_pos, pr_pos_se,
         pr_top5, pr_top5_se,
         pr_top1, pr_top1_se,
         e_y_pos, e_y_pos_se,
         e_y) %>%
  arrange(insurance_status, n_chronic_conditions_cat)

write.csv(out, "output/tables/table_twopart_profiles.csv", row.names = FALSE)
cat("CSV saved: output/tables/table_twopart_profiles.csv\n")

# -----------------------------------------------------------------------------
# 6. Formatted HTML table
# -----------------------------------------------------------------------------
cat("Rendering HTML table...\n")

fmt_pct       <- function(x, se) sprintf("%.1f%% (%.1f)", 100 * x, 100 * se)
fmt_dol       <- function(x, se) sprintf("$%s (%s)",
                                         formatC(x,  format = "d", big.mark = ","),
                                         formatC(se, format = "d", big.mark = ","))
fmt_dol_point <- function(x) sprintf("$%s",
                                     formatC(x, format = "d", big.mark = ","))

tab <- out %>%
  mutate(
    `Pr(Y > 0)`    = fmt_pct(pr_pos,  pr_pos_se),
    `Pr(Top 5%)`   = fmt_pct(pr_top5, pr_top5_se),
    `Pr(Top 1%)`   = fmt_pct(pr_top1, pr_top1_se),
    `E[Y | Y > 0]` = fmt_dol(e_y_pos, e_y_pos_se),
    `E[Y]`         = fmt_dol_point(e_y),
    Insurance      = as.character(insurance_status),
    Morbidity      = as.character(n_chronic_conditions_cat)
  ) %>%
  select(Insurance, Morbidity,
         `Pr(Y > 0)`, `Pr(Top 5%)`, `Pr(Top 1%)`,
         `E[Y | Y > 0]`, `E[Y]`)

gt_tab <- gt(tab) %>%
  tab_header(
    title    = "Two-part expected healthcare spending by profile",
    subtitle = paste0("Reference profile: Female, 55-64, White Non-Hispanic, ",
                      "College+, Middle Income, No SPD, No cancer history")
  ) %>%
  tab_source_note(
    source_note = paste0(
      "Point estimates from MEPS survey-weighted models. Standard errors in ",
      "parentheses (delta method). Uncertainty for E[Y] is not propagated ",
      "because it is constructed as the product of two model-based ",
      "predictions from separate stages of the two-part model."
    )
  ) %>%
  cols_align(align = "left",   columns = c(Insurance, Morbidity)) %>%
  cols_align(align = "right",
             columns = c(`Pr(Y > 0)`, `Pr(Top 5%)`, `Pr(Top 1%)`,
                         `E[Y | Y > 0]`, `E[Y]`))

gtsave(gt_tab, "output/tables/table_twopart_profiles.html")
cat("HTML saved: output/tables/table_twopart_profiles.html\n")

# -----------------------------------------------------------------------------
# 7. Console summary
# -----------------------------------------------------------------------------
cat("\n=============================================================\n")
cat("Two-part expected spending - 6 profiles\n")
cat("=============================================================\n")
print(
  out %>%
    transmute(
      Insurance  = insurance_status,
      Morbidity  = n_chronic_conditions_cat,
      `Pr(Y>0)`  = sprintf("%.3f", pr_pos),
      `Pr(T5%)`  = sprintf("%.3f", pr_top5),
      `Pr(T1%)`  = sprintf("%.3f", pr_top1),
      `E[Y|Y>0]` = sprintf("$%s", formatC(e_y_pos,
                                          format = "d", big.mark = ",")),
      `E[Y]`     = sprintf("$%s", formatC(e_y,
                                          format = "d", big.mark = ","))
    )
)

cat("\nStep 6 complete.\n")
