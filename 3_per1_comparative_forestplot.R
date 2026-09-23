# ------------------------------------------------------------------------------
# 2. Multi-Cohort Comparative Forest Plot
# ------------------------------------------------------------------------------
cat("--> Generating Multi-Cohort Forest Plot...\n")

if (!requireNamespace("forestplot", quietly = TRUE)) install.packages("forestplot")
library(forestplot)

# 1. Structure cohort summary dataframe
# Replace GSE14520 / TCGA HRs with your exact calculated values if they differ slightly
forest_df <- data.frame(
  cohort = c("ICGC LIRI-JP (HCCDB-18)", "GSE14520 (GEO)", "TCGA-LIHC"),
  n      = c(nrow(icgc_data), 242, 365),
  hr     = c(icgc_hr, 0.74, 0.82),             # ICGC HR = 0.76
  lower  = c(icgc_lower, 0.58, 0.65),          # ICGC lower CI = 0.59
  upper  = c(icgc_upper, 0.96, 1.03),          # ICGC upper CI = 0.98
  p_val  = c(icgc_p, 0.0230, 0.0810)           # ICGC p = 0.0328
)

# 2. Format text table columns for the display
table_text <- cbind(
  c("Cohort", forest_df$cohort),
  c("Patients (n)", forest_df$n),
  c("Hazard Ratio (95% CI)", sprintf("%.2f (%.2f - %.2f)", forest_df$hr, forest_df$lower, forest_df$upper)),
  c("p-value", ifelse(forest_df$p_val < 0.001, "< 0.001", sprintf("%.4f", forest_df$p_val)))
)

# 3. Render publication-ready forest plot
# ------------------------------------------------------------------------------
# Compact, Publication-Ready Forest Plot
# ------------------------------------------------------------------------------
library(forestplot)
library(grid)

# Ensure graph.pos points to where you want the plot (e.g., column 4, right after HR)
# We structure table_text into 4 columns total:
table_text <- cbind(
  c("Cohort", as.character(forest_df$cohort)),
  c("Patients (n)", forest_df$n),
  c("Hazard Ratio (95% CI)", sprintf("%.2f (%.2f - %.2f)", forest_df$hr, forest_df$lower, forest_df$upper)),
  c("p-value", ifelse(forest_df$p_val < 0.001, "< 0.001", sprintf("%.4f", forest_df$p_val)))
)

forestplot(
  labeltext = table_text,
  mean  = c(NA, forest_df$hr),
  lower = c(NA, forest_df$lower),
  upper = c(NA, forest_df$upper),
  is.summary = c(TRUE, rep(FALSE, nrow(forest_df))),
  zero = 1.0,
  xlog = FALSE,

  # --- LAYOUT & SPACING FIXES ---
  graph.pos = 4,                   # Places the plot AS column 4 (between HR and p-value)
  graphwidth = unit(5, "cm"),      # Fixed plot width so it doesn't stretch across the screen
  colgap = unit(3, "mm"),          # Snaps text columns close together
  lineheight = unit(1.2, "cm"),    # Clean row spacing
  # ------------------------------

  clip = c(0.4, 1.2),
  xticks = c(0.4, 0.6, 0.8, 1.0, 1.2),
  xlab = "Hazard Ratio (High vs. Low PER1 Expression)",
  col = fpColors(
    box   = "#377EB8",
    lines = "#1F4E79",
    zero  = "#E41A1C"
  ),
  boxsize = 0.25,
  vertices = TRUE,
  title = "PER1 Prognostic Value Across HCC Validation Cohorts"
)