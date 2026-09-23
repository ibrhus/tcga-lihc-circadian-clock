library(tidyverse)
library(survival)
library(survminer)

# -----------------------------------------------------------------------------
# 1. Target Genes & Expression Data Preparation
# -----------------------------------------------------------------------------
# Define target genes of interest
target_genes <- c("PER1", "CRY2", "RORC", "RORA")

# Identify available target genes in expression matrix rows
available_genes <- intersect(rownames(circadian_expr), target_genes)
message("Target genes found in dataset: ", paste(available_genes, collapse = ", "))

# Transpose matrix for target genes (Samples as rows, Genes as columns)
expr_target_df <- as.data.frame(t(circadian_expr[available_genes, , drop = FALSE])) %>%
  rownames_to_column(var = "Sample_ID")

# -----------------------------------------------------------------------------
# 2. Clean Clinical Metadata & Merge Datasets
# -----------------------------------------------------------------------------
cox_df_full <- clinical_clean %>%
  inner_join(expr_target_df, by = "Sample_ID") %>%
  mutate(
    # Clean Age
    age = as.numeric(years_to_birth),
    # Clean Gender
    gender = factor(gender),
    # Simplify Pathologic Stage (Stage I, II, III, IV)
    stage = case_when(
      grepl("stage i$", pathologic_stage, ignore.case = TRUE) ~ "Stage I",
      grepl("stage ii$", pathologic_stage, ignore.case = TRUE) ~ "Stage II",
      grepl("stage iii", pathologic_stage, ignore.case = TRUE) ~ "Stage III",
      grepl("stage iv", pathologic_stage, ignore.case = TRUE) ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    stage = factor(stage)
  )

# -----------------------------------------------------------------------------
# 3. Fit Multivariate Cox Model (Target Genes + Clinical Covariates)
# -----------------------------------------------------------------------------
covariates <- c(available_genes, "age", "gender", "stage")
multi_formula_full <- as.formula(paste("Surv(overall_survival_days, vital_status) ~", paste(covariates, collapse = " + ")))

multi_fit_clinical <- coxph(multi_formula_full, data = cox_df_full)

# Print Summary to Console
print("--- Complete Multivariate Cox Model (Target Genes + Clinical) ---")
summary(multi_fit_clinical)

# -----------------------------------------------------------------------------
# 4. Proportional Hazards Assumption Check (Diagnostics)
# -----------------------------------------------------------------------------
ph_test_clinical <- cox.zph(multi_fit_clinical)
print("--- Proportional Hazards Assumption Test (cox.zph) ---")
print(ph_test_clinical)

# -----------------------------------------------------------------------------
# 5. Export Forest Plot and Numerical Summary
# -----------------------------------------------------------------------------
# Ensure directory exists
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Generate Forest Plot
forest_p_clinical <- ggforest(
  multi_fit_clinical,
  data = cox_df_full,
  main = "Multivariate Cox Model: Target Genes & Clinical Covariates",
  fontsize = 0.8
)

# Save Forest Plot to PNG
png("results/figures/target_genes_clinical_cox_forest.png", width = 9, height = 7, units = "in", res = 300)
print(forest_p_clinical)
dev.off()

# Save Clean Numerical Summary Table to CSV
multi_summary <- as.data.frame(summary(multi_fit_clinical)$coefficients)
multi_conf <- as.data.frame(summary(multi_fit_clinical)$conf.int)

final_cox_clinical <- data.frame(
  Variable = rownames(multi_summary),
  HR = round(multi_conf$`exp(coef)`, 2),
  CI_Lower = round(multi_conf$`lower .95`, 2),
  CI_Upper = round(multi_conf$`upper .95`, 2),
  p_value = round(multi_summary$`Pr(>|z|)`, 4)
)

write.csv(final_cox_clinical, "results/multivariate_cox_target_genes_clinical.csv", row.names = FALSE)
message("Phase 5 script completed successfully. Results and plots saved to 'results/'.")