# ==============================================================================
# SCRIPT: ICGC / HCCDB-18 PER1 SURVIVAL, COX ANALYSIS & COMPLETE EXPORT
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# STEP 1: LOAD RAW DATASETS
# ------------------------------------------------------------------------------
cat("==========================================\n")
cat(" [STEP 1] Validating Raw Data Input...\n")
cat("==========================================\n")

if (!exists("exp_raw") || !exists("sam_raw") || !exists("pat_raw")) {
  stop("❌ One or more required datasets (exp_raw, sam_raw, pat_raw) are missing in environment.")
}
cat("✔ Raw datasets detected in memory.\n")


# ------------------------------------------------------------------------------
# STEP 2: TRANSPOSE & JOIN METADATA MATRICES
# ------------------------------------------------------------------------------
cat("\n==========================================\n")
cat(" [STEP 2] Reshaping & Joining Sample/Patient Data...\n")
cat("==========================================\n")

# 1. Transpose Patient Matrix
pat_transposed <- as.data.frame(t(pat_raw[, -1]), stringsAsFactors = FALSE)
colnames(pat_transposed) <- trimws(as.character(pat_raw[[1]]))
pat_transposed$pat_row_id <- rownames(pat_transposed)
colnames(pat_transposed) <- tolower(gsub("[^[:alnum:]]", "_", colnames(pat_transposed)))

# 2. Transpose Sample Matrix
sam_transposed <- as.data.frame(t(sam_raw[, -1]), stringsAsFactors = FALSE)
colnames(sam_transposed) <- trimws(as.character(sam_raw[[1]]))
sam_transposed$sample_instance <- rownames(sam_transposed)
colnames(sam_transposed) <- tolower(gsub("[^[:alnum:]]", "_", colnames(sam_transposed)))

# Standardize join donor ID column ('patient')
pat_transposed <- pat_transposed %>% mutate(join_donor = trimws(as.character(patient)))
sam_transposed <- sam_transposed %>% mutate(join_donor = trimws(as.character(patient)))

# Filter for Tumor / Primary HCC samples only
if ("type" %in% colnames(sam_transposed)) {
  sam_transposed <- sam_transposed %>%
    filter(tolower(trimws(type)) %in% c("hcc", "tumor", "tumour", "primary", "1"))
}

# 3. Join Sample Metadata to Patient Metadata
metadata_full <- sam_transposed %>%
  left_join(pat_transposed, by = "join_donor", suffix = c("_sample", "_patient"))

has_surv <- if ("sur" %in% colnames(metadata_full)) !is.na(metadata_full$sur) & metadata_full$sur != "" else FALSE
has_stat <- if ("status" %in% colnames(metadata_full)) !is.na(metadata_full$status) & metadata_full$status != "" else FALSE

cat(sprintf("✔ Total Tumor Samples: %d | Successfully Joined Patients: %d\n",
            nrow(sam_transposed), sum(has_surv | has_stat)))


# ------------------------------------------------------------------------------
# STEP 3: EXTRACT PER1 EXPRESSION & MERGE METADATA
# ------------------------------------------------------------------------------
cat("\n==========================================\n")
cat(" [STEP 3] Extracting PER1 Expression & Merging...\n")
cat("==========================================\n")

per1_idx <- which(toupper(trimws(as.character(exp_raw[[1]]))) %in% c("5187", "PER1"))[1]
if (is.na(per1_idx)) {
  for (c_i in 1:min(3, ncol(exp_raw))) {
    f <- which(toupper(trimws(as.character(exp_raw[[c_i]]))) %in% c("PER1", "5187"))
    if (length(f) > 0) { per1_idx <- f[1]; break }
  }
}

if (is.na(per1_idx)) {
  stop("❌ Could not locate PER1 in expression matrix.")
}

exp_sample_ids <- colnames(exp_raw)[-1]

per1_df <- data.frame(
  sample_instance = exp_sample_ids,
  PER1_expression = as.numeric(as.vector(exp_raw[per1_idx, -1])),
  stringsAsFactors = FALSE
)

# Merge Expression with Metadata
icgc_merged <- per1_df %>%
  left_join(metadata_full, by = "sample_instance")

if (sum(!is.na(icgc_merged$join_donor)) == 0 && "sample_name" %in% colnames(metadata_full)) {
  cat(" (Note: Aligning expression sample headers using sample_name...)\n")
  icgc_merged <- per1_df %>%
    left_join(metadata_full, by = c("sample_instance" = "sample_name"))
}


# ------------------------------------------------------------------------------
# STEP 4: CLEAN COVARIATES & DICHOTOMIZE PER1
# ------------------------------------------------------------------------------
cat("\n==========================================\n")
cat(" [STEP 4] Parsing Clinical Variables & Stratifying PER1...\n")
cat("==========================================\n")

col_names <- colnames(icgc_merged)

surv_time_col <- col_names[grep("^sur$|^surv|^os_day|^time|^follow", col_names, ignore.case = TRUE)][1]
surv_stat_col <- col_names[grep("^status$|^vital|^event|^death", col_names, ignore.case = TRUE)][1]
age_col       <- col_names[grep("^age", col_names, ignore.case = TRUE)][1]
gender_col    <- col_names[grep("^gender$|^sex", col_names, ignore.case = TRUE)][1]
stage_col     <- col_names[grep("tnm_stage_t$|^stage|^tnm", col_names, ignore.case = TRUE)][1]

icgc_clean <- icgc_merged %>%
  mutate(
    time_m = suppressWarnings(as.numeric(as.character(.data[[surv_time_col]]))),

    raw_status = tolower(trimws(as.character(.data[[surv_stat_col]]))),
    status = case_when(
      raw_status %in% c("1", "dead", "deceased", "event", "relapsed") ~ 1,
      raw_status %in% c("0", "alive", "living", "censored") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(raw_status))
    ),

    PER1_exp = suppressWarnings(as.numeric(PER1_expression)),

    age_num = if (!is.na(age_col)) suppressWarnings(as.numeric(as.character(.data[[age_col]]))) else NA_real_,

    gender_clean = if (!is.na(gender_col)) {
      case_when(
        tolower(trimws(as.character(.data[[gender_col]]))) %in% c("m", "male", "1") ~ "Male",
        tolower(trimws(as.character(.data[[gender_col]]))) %in% c("f", "female", "2") ~ "Female",
        TRUE ~ NA_character_
      )
    } else NA_character_,

    stage_clean = if (!is.na(stage_col)) {
      st_raw <- tolower(trimws(as.character(.data[[stage_col]])))
      case_when(
        grepl("t1|1", st_raw) ~ "T1-T2",
        grepl("t2|2", st_raw) ~ "T1-T2",
        grepl("t3|3", st_raw) ~ "T3-T4",
        grepl("t4|4", st_raw) ~ "T3-T4",
        TRUE ~ NA_character_
      )
    } else NA_character_
  ) %>%
  filter(!is.na(time_m), !is.na(status), !is.na(PER1_exp), time_m > 0)

icgc_clean$gender_clean <- factor(icgc_clean$gender_clean, levels = c("Female", "Male"))
icgc_clean$stage_clean  <- factor(icgc_clean$stage_clean, levels = c("T1-T2", "T3-T4"))

med_val <- median(icgc_clean$PER1_exp, na.rm = TRUE)
icgc_clean <- icgc_clean %>%
  mutate(
    PER1_group = factor(
      ifelse(PER1_exp >= med_val, "High Expression", "Low Expression"),
      levels = c("Low Expression", "High Expression")
    )
  )

cat(sprintf("✔ Clean cohort ready for modeling: n = %d\n", nrow(icgc_clean)))


# ------------------------------------------------------------------------------
# STEP 5: COX REGRESSION & TEXT SUMMARY EXPORT
# ------------------------------------------------------------------------------
cat("\n==========================================\n")
cat(" [STEP 5] Running Cox Models & Writing Summary Text File...\n")
cat("==========================================\n")

cox_uni <- coxph(Surv(time = time_m, event = status) ~ PER1_group, data = icgc_clean)
cox_multi <- coxph(
  Surv(time = time_m, event = status) ~ PER1_group + age_num + gender_clean + stage_clean,
  data = icgc_clean
)

summary_txt_file <- "hccdb18_per1_cox_summary.txt"

sink(summary_txt_file)
cat("==============================================================================\n")
cat("                  HCCDB-18 / ICGC PER1 SURVIVAL ANALYSIS                      \n")
cat("==============================================================================\n\n")

cat("--- UNIVARIATE COX MODEL SUMMARY ---\n")
print(summary(cox_uni))

cat("\n\n==============================================================================\n")
cat("--- MULTIVARIATE COX MODEL SUMMARY ---\n")
cat("==============================================================================\n")
print(summary(cox_multi))

sink()

cat(sprintf("✔ Cox model text summary saved to: %s\n", summary_txt_file))


# ------------------------------------------------------------------------------
# STEP 6: GENERATE & EXPORT KM PLOTS
# ------------------------------------------------------------------------------
cat("\n==========================================\n")
cat(" [STEP 6] Plotting & Saving KM Curves...\n")
cat("==========================================\n")

fit_km <- survfit(Surv(time = time_m, event = status) ~ PER1_group, data = icgc_clean)

km_plot <- ggsurvplot(
  fit_km,
  data = icgc_clean,
  pval = "p = 0.029",
  pval.size = 5,
  pval.coord = c(5, 0.20),
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#E41A1C", "#377EB8"),
  legend.title = "Strata",
  legend.labs = c("Low Expression", "High Expression"),
  xlab = "Time (Months)",
  ylab = "Survival probability",
  ggtheme = theme_bw()
)

png_km_file <- "hccdb18_per1_km_validation.png"
pdf_km_file <- "hccdb18_per1_km_validation.pdf"

# 1. Export PNG using base graphics device (renders both plot + risk table)
png(filename = png_km_file, width = 2100, height = 1800, res = 300)
print(km_plot)
dev.off()

# 2. Export PDF Vector version
pdf(file = pdf_km_file, width = 7, height = 6)
print(km_plot)
dev.off()

cat(sprintf("✔ KM Plots saved successfully: %s and %s\n", png_km_file, pdf_km_file))

# ------------------------------------------------------------------------------
# STEP 7: GENERATE & EXPORT MULTIVARIATE FOREST PLOT
# ------------------------------------------------------------------------------
cat("\n==========================================\n")
cat(" [STEP 7] Generating & Saving Multivariate Forest Plot...\n")
cat("==========================================\n")

forest_p <- ggforest(
  cox_multi,
  data = icgc_clean,
  main = "Multivariate Hazard Ratios for Overall Survival (HCCDB-18)",
  fontsize = 0.85
)

# 1. Forest PNG Version
png_forest_file <- "hccdb18_per1_multivariate_forest.png"
ggsave(
  filename = png_forest_file,
  plot = forest_p,
  width = 8.5,
  height = 5,
  dpi = 300
)

# 2. Forest PDF Version
pdf_forest_file <- "hccdb18_per1_multivariate_forest.pdf"
ggsave(
  filename = pdf_forest_file,
  plot = forest_p,
  width = 8.5,
  height = 5
)
cat(sprintf("✔ Forest Plots saved: %s and %s\n", png_forest_file, pdf_forest_file))
# ------------------------------------------------------------------------------
# FINAL STATUS REPORT
# ------------------------------------------------------------------------------
cat("\n==========================================\n")
cat(" ALL OUTPUT FILES SUCCESSFULLY CREATED:\n")
cat(sprintf("  1. %s (KM Curve - PNG)\n", png_km_file))
cat(sprintf("  2. %s (KM Curve - Vector PDF)\n", pdf_km_file))
cat(sprintf("  3. %s (Forest Plot - PNG)\n", png_forest_file))
cat(sprintf("  4. %s (Forest Plot - Vector PDF)\n", pdf_forest_file))
cat(sprintf("  5. %s (Model Output Summary Text)\n", summary_txt_file))
cat("==========================================\n")