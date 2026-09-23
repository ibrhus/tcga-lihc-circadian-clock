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
# 1. Download & Load GSE14520 Expression Matrix
# ==============================================================================
cat("--> [1/6] Loading GSE14520 Series Matrix from GEO...\n")

gse <- getGEO("GSE14520", destdir = data_dir, GSEMatrix = TRUE)

# Select GPL3921 platform (HT Human Genome U133A Array with primary survival cohort)
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

# Identify dynamic column names
t_col <- grep("survival.*month|os.*time|survival_months|survival.months", colnames(clin_raw), value = TRUE)[1]
s_col <- grep("survival.*status|os.*status|survival_status|survival.status", colnames(clin_raw), value = TRUE)[1]

cat("--> Using Time Col:", t_col, "| Status Col:", s_col, "\n")

clin_clean <- clin_raw %>%
  filter(!is.na(.data[[t_col]]), !is.na(.data[[s_col]])) %>%
  mutate(
    time_m   = as.numeric(as.character(.data[[t_col]])),
    raw_stat = as.numeric(as.character(.data[[s_col]])),
    time     = time_m * 30.4375, # Convert months to days
    status   = ifelse(raw_stat %in% c(1, "1", "Dead", "dead"), 1, 0),

    # Clean keys for matching: remove hyphens, underscores, and spaces (e.g. LCS_193A -> LCS193A)
    gsm_key  = str_trim(as.character(affy_gsm)),
    lcs_key  = gsub("[^A-Za-z0-9]", "", as.character(lcs.id)),

    age_clean = if ("age" %in% colnames(.)) as.numeric(as.character(age)) else NA_real_,
    stage_raw = if ("tnm.stage" %in% colnames(.)) as.character(tnm.stage) else if ("stage" %in% colnames(.)) as.character(stage) else NA_character_
  )

# Filter for Tumor samples only
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
    # Standardize sample title (e.g. "Liver Tumor Tissue LCS-079A" -> "LCS079A")
    lcs_key = gsub("[^A-Za-z0-9]", "", str_extract(title, "LCS[-_ ]?[0-9]+[A-Za-z]?"))
  ) %>%
  # Filter out non-tumor samples from expression set
  filter(grepl("Tumor", title, ignore.case = TRUE) & !grepl("Non-Tumor", title, ignore.case = TRUE))

# Primary Strategy: Join on GSM ID (geo_accession <-> affy_gsm)
geo_data <- expr_df %>% inner_join(clin_clean, by = "gsm_key")

# Secondary Strategy: Join on Normalized LCS ID
if (nrow(geo_data) == 0) {
  cat("--> GSM join returned 0 rows. Falling back to normalized LCS sample ID join...\n")
  geo_data <- expr_df %>% inner_join(clin_clean, by = "lcs_key")
}

# Filter out missing or non-positive survival times
geo_data <- geo_data %>%
  filter(!is.na(time), !is.na(status), !is.na(PER1_expression), time > 0)

cat("--> Final Cohort Check:\n")
cat("    - Valid Tumor Patients:", nrow(geo_data), "\n")
cat("    - Death Events         :", sum(geo_data$status == 1), "\n")

if (nrow(geo_data) == 0) {
  stop("FATAL: Could not match sample IDs between expression matrix and clinical table.")
}

# ==============================================================================
# 5. Kaplan-Meier Survival Analysis
# ==============================================================================
cat("--> [5/6] Fitting Kaplan-Meier Survival Analysis...\n")

geo_data <- geo_data %>%
  mutate(
    PER1_group = ifelse(PER1_expression >= median(PER1_expression, na.rm = TRUE),
                        "High Expression", "Low Expression")
  )

km_fit <- survfit(Surv(time, status) ~ PER1_group, data = geo_data)

p_km <- ggsurvplot(
  km_fit,
  data = geo_data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.title = "PER1 Expression",
  legend.labs = c("High Expression", "Low Expression"),
  palette = c("#E41A1C", "#377EB8"),
  title = "External Validation: PER1 Overall Survival (GSE14520)",
  xlab = "Time (Days)",
  ylab = "Overall Survival Probability",
  ggtheme = theme_bw(base_size = 12)
)

print(p_km)

ggsave(
  filename = file.path(results_dir, "gse14520_per1_km_validation.png"),
  plot = p_km$plot,
  width = 7,
  height = 6,
  dpi = 300
)
# ==============================================================================
# 6. Multivariate Cox Regression Analysis
# ==============================================================================
cat("--> [6/6] Fitting Multivariate Cox Model...\n")

covars <- c("PER1_expression")
if (sum(!is.na(geo_data$age_clean)) > 10) covars <- c(covars, "age_clean")
if ("stage_raw" %in% colnames(geo_data) && length(unique(na.omit(geo_data$stage_raw))) > 1) {
  covars <- c(covars, "stage_raw")
}

cox_formula <- as.formula(paste("Surv(time, status) ~", paste(covars, collapse = " + ")))
cox_gse     <- coxph(cox_formula, data = geo_data)

cat("\n=================== GSE14520 MULTIVARIATE COX MODEL ===================\n")
print(summary(cox_gse))
cat("========================================================================\n\n")

# Corrected ggforest call (removed unused censor = TRUE argument)
p_forest <- ggforest(
  cox_gse,
  data = geo_data,
  main = "Multivariate Cox Model: GSE14520 External Validation"
)

ggsave(
  filename = file.path(results_dir, "gse14520_per1_multivariate_forest.png"),
  plot = p_forest,
  width = 8,
  height = 5,
  dpi = 300
)

cat("--> Pipeline finished successfully! Figures saved to 'results/'.\n")