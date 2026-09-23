# ==============================================================================
# Project: TCGA-LIHC Circadian Clock Analysis
# Phase 3: Differential Expression Analysis (DEA) of Circadian Panel
# Author: Hussam Ibrahim
# Date: July 2026
# ==============================================================================

# --- 1. Load Required Libraries ---
library(curatedTCGAData)
library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(ggrepel)

# --- 2. Create Directory Structure ---
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# --- 3. Define Circadian Target Panel (21 Genes) ---
target_genes <- c(
  "ARNTL", "BHLHE40", "BHLHE41", "CDKN1A", "CLOCK", "CRY1", "CRY2",
  "HIF1A", "HNF4A", "NPAS2", "NR1D1", "NR1D2", "NR1H4", "NR1I3",
  "PER1", "PER2", "PER3", "RORA", "RORB", "RORC", "WEE1"
)

# --- 4. Load TCGA-LIHC RNA-Seq Data ---
message("Fetching TCGA-LIHC RNA-Seq dataset...")
lihc_rna <- curatedTCGAData(
  diseaseCode = "LIHC",
  assays = "RNASeq2GeneNorm",
  version = "2.1.1",
  dry.run = FALSE
)

se_rna <- lihc_rna[[1]]

# --- 5. Clean Samples and Construct Metadata ---
# Extract sample codes: '01' = Primary Solid Tumor, '11' = Solid Tissue Normal
sample_ids <- colnames(se_rna)
valid_samples <- sample_ids[substr(sample_ids, 14, 15) %in% c("01", "11")]

# Extract matrix for valid samples
expr_matrix <- assays(se_rna)[[1]][, valid_samples]

# Build sample metadata table
col_data <- data.frame(
  Sample_ID = valid_samples,
  Condition = ifelse(substr(valid_samples, 14, 15) == "01", "Tumor", "Normal"),
  row.names = valid_samples
)
col_data$Condition <- factor(col_data$Condition, levels = c("Normal", "Tumor"))

# Filter expression matrix for the 21 target circadian genes
circadian_expr <- expr_matrix[rownames(expr_matrix) %in% target_genes, , drop = FALSE]

message("Cohort summary: ", sum(col_data$Condition == "Tumor"), " Primary Tumors, ",
        sum(col_data$Condition == "Normal"), " Adjacent Normals across ",
        nrow(circadian_expr), " panel genes.")

# --- 6. Differential Expression Analysis with limma ---
design <- model.matrix(~ Condition, data = col_data)
fit <- lmFit(circadian_expr, design)
fit <- eBayes(fit)

# Extract and annotate results
dea_results <- topTable(fit, coef = "ConditionTumor", number = Inf) %>%
  rownames_to_column(var = "Hugo_Symbol") %>%
  as_tibble()

dea_summary <- dea_results %>%
  mutate(Status = case_when(
    adj.P.Val < 0.05 & logFC > 0.5 ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < -0.5 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Save statistical summary table
write_csv(dea_summary, "data/processed/circadian_dea_results.csv")
message("Saved DEA results to 'data/processed/circadian_dea_results.csv'")

# --- 7. Plot 1: Expression Boxplots (Top Downregulated Genes) ---
top_genes <- c("PER1", "PER2", "BHLHE40", "NR1I3", "RORA", "WEE1")

plot_df <- circadian_expr[top_genes, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Hugo_Symbol") %>%
  pivot_longer(-Hugo_Symbol, names_to = "Sample_ID", values_to = "Expression") %>%
  left_join(col_data, by = "Sample_ID")

p_box <- ggplot(plot_df, aes(x = Hugo_Symbol, y = Expression, fill = Condition)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("Normal" = "#2b8cbe", "Tumor" = "#de2d26")) +
  theme_bw() +
  labs(
    title = "Circadian Gene Expression: Tumor vs. Adjacent Normal (TCGA-LIHC)",
    x = "Gene Symbol",
    y = "Normalized Expression (RSEM)",
    fill = "Condition"
  ) +
  theme(text = element_text(size = 12))

png("results/figures/circadian_dea_boxplots.png", width = 10, height = 6, units = "in", res = 300)
print(p_box)
dev.off()
message("Saved boxplots to 'results/figures/circadian_dea_boxplots.png'")

# --- 8. Plot 2: Volcano Plot (21-Gene Panel) ---
volcano_df <- dea_summary %>%
  mutate(
    Significance = case_when(
      adj.P.Val < 0.05 & logFC > 0.5 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Hugo_Symbol), size = 3.5, max.overlaps = 20) +
  scale_color_manual(values = c("Downregulated" = "#2b8cbe", "Upregulated" = "#de2d26", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  theme_bw() +
  labs(
    title = "Volcano Plot: Circadian Panel DEA (TCGA-LIHC)",
    x = "Log2 Fold Change (Tumor vs. Normal)",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme(text = element_text(size = 12))

png("results/figures/circadian_dea_volcano.png", width = 8, height = 6, units = "in", res = 300)
print(p_volcano)
dev.off()
message("Saved volcano plot to 'results/figures/circadian_dea_volcano.png'")

# --- Complete ---
message("Phase 3 Script Execution Complete!")