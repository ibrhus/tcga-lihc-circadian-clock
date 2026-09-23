# ==============================================================================
# Project: TCGA-LIHC Circadian Clock Analysis
# Phase 4: Survival & Clinical Correlation Analysis
# Author: Hussam Ibrahim
# Date: July 2026
# ==============================================================================

# --- 1. Load Required Libraries ---
library(curatedTCGAData)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)

# --- 2. Create Output Directory Structure ---
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/km_plots", recursive = TRUE, showWarnings = FALSE)

# --- 3. Extract & Clean Clinical Metadata ---
message("Processing TCGA-LIHC clinical survival metadata...")
clinical_raw <- as.data.frame(colData(lihc_rna)) %>%
  rownames_to_column(var = "Patient_ID")

clinical_clean <- clinical_raw %>%
  mutate(
    # Vital status: 1 = Deceased, 0 = Living/Censored
    status = ifelse(tolower(as.character(vital_status)) %in% c("dead", "1"), 1, 0),
    days_death = as.numeric(days_to_death),
    days_follow = as.numeric(days_to_last_followup),
    # Compute overall survival time in days
    time = ifelse(!is.na(days_death), days_death, days_follow)
  ) %>%
  filter(!is.na(time) & time > 0)

# --- 4. Link Primary Tumor RNA-Seq Cohort to Survival Data ---
# Extract Primary Solid Tumor barcodes ('01')
tumor_sample_ids <- colnames(circadian_expr)[substr(colnames(circadian_expr), 14, 15) == "01"]
tumor_expr <- circadian_expr[, tumor_sample_ids, drop = FALSE]

# Map 12-character Patient IDs to 16+ character Sample Barcodes
sample_patient_map <- data.frame(
  Sample_ID = tumor_sample_ids,
  Patient_ID = substr(tumor_sample_ids, 1, 12),
  stringsAsFactors = FALSE
)

# Merge expression and clinical records (keep distinct patient records)
surv_df <- sample_patient_map %>%
  inner_join(clinical_clean, by = "Patient_ID") %>%
  distinct(Patient_ID, .keep_all = TRUE)

message("Cohort linked: ", nrow(surv_df), " primary tumor samples matched with valid survival outcomes.")

# --- 5. Batch Cox Proportional Hazards & Kaplan-Meier Analysis ---
target_surv_genes <- c("PER1", "PER2", "BHLHE40", "NR1I3", "RORA", "WEE1", "RORC", "CRY1", "CRY2")

cox_results <- list()

for (gene in target_surv_genes) {
  group_col <- paste0(gene, "_group")

  # Stratify patients by median expression level into High vs. Low cohorts
  g_vals <- tumor_expr[gene, surv_df$Sample_ID]
  med_val <- median(g_vals, na.rm = TRUE)
  surv_df[[group_col]] <- ifelse(g_vals >= med_val, "High", "Low")

  # Construct dynamic formulas
  f_str <- paste0("Surv(time, status) ~ ", group_col)
  f_form <- as.formula(f_str)

  # Fit Cox Proportional Hazards Regression
  cox_fit <- coxph(f_form, data = surv_df)
  cox_sum <- summary(cox_fit)

  cox_results[[gene]] <- tibble(
    Hugo_Symbol = gene,
    Hazard_Ratio = cox_sum$coefficients[1, "exp(coef)"],
    HR_Lower_95 = cox_sum$conf.int[1, "lower .95"],
    HR_Upper_95 = cox_sum$conf.int[1, "upper .95"],
    P_Value = cox_sum$coefficients[1, "Pr(>|z|)"]
  )

  # Fit Kaplan-Meier Curve (eval/bquote locks environment for ggsurvplot)
  fit_km <- eval(bquote(survfit(.(f_form), data = surv_df)))

  # Generate KM Survival Plot
  km_plot <- ggsurvplot(
    fit_km,
    data = surv_df,
    pval = TRUE,
    pval.size = 5,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    ggtheme = theme_bw(),
    palette = c("#de2d26", "#2b8cbe"),
    title = paste0("Overall Survival: ", gene, " Expression in TCGA-LIHC"),
    xlab = "Time (Days)",
    legend.labs = c("High Expression", "Low Expression")
  )

  # Save PNG Figure
  png(paste0("results/figures/km_plots/", gene, "_survival_km.png"),
      width = 8, height = 7, units = "in", res = 300)
  print(km_plot)
  dev.off()
}

# --- 6. Export Statistical Summary ---
cox_summary_df <- bind_rows(cox_results)
write_csv(cox_summary_df, "data/processed/circadian_survival_cox_summary.csv")

message("Saved Cox summary to 'data/processed/circadian_survival_cox_summary.csv'")
message("Phase 4 Execution Complete! All KM plots exported to 'results/figures/km_plots/'.")