# ------------------------------------------------------------------------------
# 1. Multivariate Cox Proportional Hazards Model
# ------------------------------------------------------------------------------
cat("--> [1/2] Preparing clinical covariates & running Cox regression...\n")

# Clean clinical covariates using your exact column names
icgc_data <- icgc_data %>%
  mutate(
    age_num      = suppressWarnings(as.numeric(as.character(age))),
    gender_clean = tolower(trimws(as.character(gender))),
    stage_clean  = tolower(trimws(as.character(tnm_stage_t)))
  )

# Fit Multivariate Cox Model
cox_multi <- coxph(
  Surv(time = time_m, event = status) ~ PER1_expression + age_num + gender_clean + stage_clean,
  data = icgc_data
)

# Display complete summary
summary(cox_multi)

# Extract ICGC PER1 Adjusted Hazard Ratio & 95% CI
cox_icgc_sum <- summary(cox_multi)
icgc_hr    <- cox_icgc_sum$coefficients["PER1_expression", "exp(coef)"]
icgc_lower <- cox_icgc_sum$conf.int["PER1_expression", "lower .95"]
icgc_upper <- cox_icgc_sum$conf.int["PER1_expression", "upper .95"]
icgc_p     <- cox_icgc_sum$coefficients["PER1_expression", "Pr(>|z|)"]

cat("\n------------------------------------------------------------------\n")
cat(sprintf("ICGC PER1 Adjusted HR: %.2f (95%% CI: %.2f - %.2f), p = %.4f\n",
            icgc_hr, icgc_lower, icgc_upper, icgc_p))
cat("------------------------------------------------------------------\n")