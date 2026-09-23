# ==============================================================================
# Script: PER1 Survival Analysis (HCCDB-18 / ICGC LIRI-JP)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
})

cat("--> [1/5] Libraries loaded successfully.\n")

# ------------------------------------------------------------------------------
# 2. Transpose Metadata Tables & Standardize Key Names
# ------------------------------------------------------------------------------
cat("--> [2/5] Transposing clinical tables...\n")

# Transpose patient_df (Row headers -> patient_id column)
patient_clean <- patient_df %>%
  tibble::column_to_rownames(var = colnames(patient_df)[1]) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "patient_id")

colnames(patient_clean) <- tolower(trimws(colnames(patient_clean)))
patient_clean$patient_id <- tolower(trimws(patient_clean$patient_id))

# Transpose sample_df (Row headers -> sample_col_id column)
sample_clean <- sample_df %>%
  tibble::column_to_rownames(var = colnames(sample_df)[1]) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample_col_id")

colnames(sample_clean) <- tolower(trimws(colnames(sample_clean)))
sample_clean$sample_col_id <- tolower(trimws(sample_clean$sample_col_id))
if ("patient_id" %in% colnames(sample_clean)) {
  sample_clean$patient_id <- tolower(trimws(sample_clean$patient_id))
}

# Standardize PER1 expression table keys
per1_clean <- per1_df %>%
  mutate(sample_id = tolower(trimws(as.character(sample_id))))

# ------------------------------------------------------------------------------
# 3. Merge Clinical Metadata with PER1 Expression Data & Clean Variables
# ------------------------------------------------------------------------------
cat("--> [3/5] Merging metadata and cleaning survival metrics...\n")

# Merge sample and patient metadata on patient_id
if ("patient_id" %in% colnames(sample_clean) && "patient_id" %in% colnames(patient_clean)) {
  clin_merged <- inner_join(sample_clean, patient_clean, by = "patient_id")
} else {
  clin_merged <- bind_cols(sample_clean, patient_clean)
}

# Merge with per1_clean using sample_col_id <-> sample_id
icgc_data <- clin_merged %>%
  inner_join(per1_clean, by = c("sample_col_id" = "sample_id")) %>%
  mutate(
    # Clean time (sur is already in months)
    raw_time = suppressWarnings(as.numeric(as.character(sur))),

    # Map status safely across numeric/string formats
    status_clean = tolower(trimws(as.character(status))),
    status = case_when(
      status_clean %in% c("1", "dead", "deceased", "event") ~ 1,
      status_clean %in% c("0", "alive", "living", "censored") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(status_clean))
    ),

    PER1_expression = suppressWarnings(as.numeric(as.character(PER1_expression)))
  ) %>%
  filter(!is.na(raw_time), !is.na(status), !is.na(PER1_expression), raw_time > 0) %>%
  mutate(
    time_m = raw_time,
    per1_group = ifelse(PER1_expression >= median(PER1_expression, na.rm = TRUE), "High", "Low")
  )

cat("--> Final Cohort Summary:\n")
cat("    - Valid Patients:", nrow(icgc_data), "\n")
cat("    - Death Events  :", sum(icgc_data$status == 1, na.rm = TRUE), "\n")

# ------------------------------------------------------------------------------
# 4. Survival Modeling (Kaplan-Meier & Cox Proportional Hazards)
# ------------------------------------------------------------------------------
cat("--> [4/5] Fitting survival models...\n")

# Fit Kaplan-Meier curve
surv_obj <- Surv(time = icgc_data$time_m, event = icgc_data$status)
fit_km   <- survfit(surv_obj ~ per1_group, data = icgc_data)

# Fit Cox regression model
cox_fit  <- coxph(surv_obj ~ PER1_expression, data = icgc_data)
cox_sum  <- summary(cox_fit)

hr     <- round(cox_sum$coefficients[1, "exp(coef)"], 2)
p_val  <- format.pval(cox_sum$coefficients[1, "Pr(>|z|)"], digits = 3)

cat("--> Cox Model Results:\n")
cat("    - Hazard Ratio (HR):", hr, "\n")
cat("    - P-value          :", p_val, "\n")

# ------------------------------------------------------------------------------
# 5. Generate & Display Kaplan-Meier Plot
# ------------------------------------------------------------------------------
cat("--> [5/5] Generating KM Plot...\n")

km_plot <- ggsurvplot(
  fit_km,
  data = icgc_data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#E41A1C", "#377EB8"),
  legend.labs = c("PER1 High", "PER1 Low"),
  legend.title = "", # Left empty to avoid ggplot2 scale title warnings
  xlab = "Time (Months)",
  ylab = "Overall Survival Probability",
  title = "Overall Survival by PER1 Expression (HCCDB-18)",
  ggtheme = theme_classic()
)

# Explicitly apply legend titles to ggplot component
km_plot$plot <- km_plot$plot +
  guides(color = guide_legend(title = "Expression"),
         fill = guide_legend(title = "Expression"))

print(km_plot)