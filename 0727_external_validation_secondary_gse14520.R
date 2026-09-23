# ==============================================================================
# Script: 02_external_validation_per1.R
# Aim: External Validation of PER1 Prognostic Value in GSE14520
# Author: Hussam Ibrahim
# ==============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(tidyverse)
  library(GEOquery)
})

# Define Directories
data_dir    <- "data_external/gse14520"
results_dir <- "results"

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# ==============================================================================
# 1. Download & Load GSE14520 Series Matrix from GEO
# ==============================================================================
cat("--> [1/6] Loading GSE14520 Series Matrix from GEO...\n")

gse <- getGEO("GSE14520", destdir = data_dir, GSEMatrix = TRUE)

if (length(gse) > 1) {
  gpl_names  <- sapply(gse, function(x) annotation(x))
  target_idx <- grep("GPL3921", gpl_names)
  if (length(target_idx) == 0) target_idx <- 1
  gse_set    <- gse[[target_idx[1]]]
} else {
  gse_set <- gse[[1]]
}

pdata_df   <- pData(gse_set)
expr_mat   <- exprs(gse_set)
feature_df <- fData(gse_set)

cat("--> Loaded Expression Matrix for", ncol(expr_mat), "samples.\n")

# ==============================================================================
# 2. Extract PER1 Probe Expression
# ==============================================================================
cat("--> [2/6] Extracting PER1 probe expression...\n")

per1_probes <- feature_df %>%
  filter(grepl("^PER1$", `Gene Symbol`) | `Gene Symbol` == "PER1") %>%
  rownames()

if (length(per1_probes) == 0) {
  per1_probes <- rownames(feature_df)[grep("PER1", feature_df$`Gene Symbol`, ignore.case = TRUE)]
}

if (length(per1_probes) == 0) {
  stop("FATAL: PER1 probe not found in feature data.")
}

if (length(per1_probes) > 1) {
  per1_expr_vec <- colMeans(expr_mat[per1_probes, , drop = FALSE], na.rm = TRUE)
} else {
  per1_expr_vec <- expr_mat[per1_probes, ]
}

# ==============================================================================
# 3. Read & Clean Clinical Metadata Table
# ==============================================================================
cat("--> [3/6] Reading clinical supplementary metadata...\n")

supp_files <- list.files(data_dir, pattern = "Extra|Clinical|Supplement", full.names = TRUE)
supp_files <- supp_files[!grepl("\\.tar$", supp_files, ignore.case = TRUE)]

if (length(supp_files) == 0) {
  cat("--> Downloading clinical metadata file from GEO...\n")
  getGEOSuppFiles("GSE14520", baseDir = data_dir, makeDirectory = FALSE, filter_regex = "Extra|Clinical|Supplement")
  supp_files <- list.files(data_dir, pattern = "Extra|Clinical|Supplement", full.names = TRUE)
  supp_files <- supp_files[!grepl("\\.tar$", supp_files, ignore.case = TRUE)]
}

clin_file <- supp_files[1]
cat("--> Parsing metadata file:", basename(clin_file), "\n")

clin_raw <- read.delim(clin_file, header = TRUE, stringsAsFactors = FALSE)
colnames(clin_raw) <- tolower(colnames(clin_raw))

t_col <- grep("survival.*month|os.*time|survival_months|survival.months", colnames(clin_raw), value = TRUE)[1]
s_col <- grep("survival.*status|os.*status|survival_status|survival.status", colnames(clin_raw), value = TRUE)[1]

cat("--> Using Time Col:", t_col, "| Status Col:", s_col, "\n")

clin_clean <- clin_raw %>%
  filter(!is.na(.data[[t_col]]), !is.na(.data[[s_col]])) %>%
  mutate(
    time_m   = as.numeric(as.character(.data[[t_col]])),
    raw_stat = as.numeric(as.character(.data[[s_col]])),
    time     = time_m * 30.4375,
    status   = ifelse(raw_stat %in% c(1, "1", "Dead", "dead"), 1, 0),

    gsm_key  = str_trim(as.character(affy_gsm)),
    lcs_key  = gsub("[^A-Za-z0-9]", "", as.character(lcs.id)),

    age_clean = if ("age" %in% colnames(.)) as.numeric(as.character(age)) else NA_real_,
    stage_raw = if ("tnm.stage" %in% colnames(.)) as.character(tnm.stage) else if ("stage" %in% colnames(.)) as.character(stage) else NA_character_
  )

if ("tissue.type" %in% colnames(clin_clean)) {
  clin_clean <- clin_clean %>% filter(tolower(tissue.type) == "tumor")
}

# ==============================================================================
# 4. Filter Expression Matrix & Merge Data
# ==============================================================================
cat("--> [4/6] Merging expression & clinical metadata...\n")

expr_df <- data.frame(
  geo_accession = rownames(pdata_df),
  title         = as.character(pdata_df$title),
  PER1_expression = as.numeric(per1_expr_vec),
  stringsAsFactors = FALSE
) %>%
  mutate(
    gsm_key = str_trim(geo_accession),
    lcs_key = gsub("[^A-Za-z0-9]", "", str_extract(title, "LCS[-_ ]?[0-9]+[A-Za-z]?"))
  ) %>%
  filter(grepl("Tumor", title, ignore.case = TRUE) & !grepl("Non-Tumor", title, ignore.case = TRUE))

geo_data <- expr_df %>% inner_join(clin_clean, by = "gsm_key")

if (nrow(geo_data) == 0) {
  cat("--> GSM join returned 0 rows. Falling back to normalized LCS sample ID join...\n")
  geo_data <- expr_df %>% inner_join(clin_clean, by = "lcs_key")
}

geo_data <- geo_data %>%
  filter(!is.na(time), !is.na(status), !is.na(PER1_expression), time > 0)

cat("--> Final Cohort Check:\n")
cat("    - Valid Tumor Patients:", nrow(geo_data), "\n")
cat("    - Death Events         :", sum(geo_data$status == 1), "\n")

if (nrow(geo_data) == 0) {
  stop("FATAL: Could not match sample IDs between expression matrix and clinical table.")
}

# ==============================================================================
# 5. Univariate & Multivariate Cox Models + Text Output Export
# ==============================================================================
cat("--> [5/6] Fitting Survival Models & Exporting Text Summary...\n")

geo_data <- geo_data %>%
  mutate(
    PER1_group = factor(
      ifelse(PER1_expression >= median(PER1_expression, na.rm = TRUE), "High Expression", "Low Expression"),
      levels = c("Low Expression", "High Expression")
    )
  )

cox_uni <- coxph(Surv(time, status) ~ PER1_group, data = geo_data)

covars <- c("PER1_group")
if (sum(!is.na(geo_data$age_clean)) > 10) covars <- c(covars, "age_clean")
if ("stage_raw" %in% colnames(geo_data) && length(unique(na.omit(geo_data$stage_raw))) > 1) {
  covars <- c(covars, "stage_raw")
}

cox_formula <- as.formula(paste("Surv(time, status) ~", paste(covars, collapse = " + ")))
cox_multi   <- coxph(cox_formula, data = geo_data)

# Export exact outputs to text summary file
summary_txt_file <- file.path(results_dir, "gse14520_per1_cox_summary.txt")
sink(summary_txt_file)
cat("==============================================================================\n")
cat("                  GSE14520 PER1 EXTERNAL VALIDATION ANALYSIS                  \n")
cat("==============================================================================\n\n")

cat("--- UNIVARIATE COX MODEL SUMMARY ---\n")
print(summary(cox_uni))

cat("\n\n==============================================================================\n")
cat("--- MULTIVARIATE COX MODEL SUMMARY ---\n")
cat("==============================================================================\n")
print(summary(cox_multi))
sink()

cat(sprintf("✔ Cox model text summary saved to: %s\n", summary_txt_file))

# ==============================================================================
# 6. Generate & Export Clean KM Plots & Forest Plots (PNG & PDF)
# ==============================================================================
cat("--> [6/6] Rendering and Saving Figures (PNG & PDF formats)...\n")

km_fit <- survfit(Surv(time, status) ~ PER1_group, data = geo_data)

p_km <- ggsurvplot(
  km_fit,
  data = geo_data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,            # Set to FALSE for clean unshaded curves
  risk.table = TRUE,
  legend.title = "PER1 Expression",
  legend.labs = c("Low Expression", "High Expression"),
  palette = c("#E41A1C", "#377EB8"),
  title = "External Validation: PER1 Overall Survival (GSE14520)",
  xlab = "Time (Days)",
  ylab = "Overall Survival Probability",
  ggtheme = theme_bw(base_size = 12)
)

# 1. Export KM PNG Version
png_km_file <- file.path(results_dir, "gse14520_per1_km_validation.png")
png(filename = png_km_file, width = 2100, height = 1800, res = 300)
print(p_km)
dev.off()

# 2. Export KM PDF Version
pdf_km_file <- file.path(results_dir, "gse14520_per1_km_validation.pdf")
pdf(file = pdf_km_file, width = 7, height = 6)
print(p_km)
dev.off()

# Multivariate Forest Plot
p_forest <- ggforest(
  cox_multi,
  data = geo_data,
  main = "Multivariate Hazard Ratios: GSE14520 External Validation",
  fontsize = 0.85
)

# 3. Export Forest PNG Version
png_forest_file <- file.path(results_dir, "gse14520_per1_multivariate_forest.png")
ggsave(filename = png_forest_file, plot = p_forest, width = 8.5, height = 5, dpi = 300)

# 4. Export Forest PDF Version
pdf_forest_file <- file.path(results_dir, "gse14520_per1_multivariate_forest.pdf")
ggsave(filename = pdf_forest_file, plot = p_forest, width = 8.5, height = 5)

cat("\n==========================================\n")
cat(" ALL GSE14520 EXTERNAL VALIDATION FILES CREATED:\n")
cat(sprintf("  1. %s (KM Curve - PNG)\n", png_km_file))
cat(sprintf("  2. %s (KM Curve - Vector PDF)\n", pdf_km_file))
cat(sprintf("  3. %s (Forest Plot - PNG)\n", png_forest_file))
cat(sprintf("  4. %s (Forest Plot - Vector PDF)\n", pdf_forest_file))
cat(sprintf("  5. %s (Model Output Summary Text)\n", summary_txt_file))
cat("==========================================\n")