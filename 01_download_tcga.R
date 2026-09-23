# --------------------------------------------------
# 01_download_and_plot_tcga.R
# --------------------------------------------------

# 1. Set working directory
setwd("/Users/husibr/hcc_circadian_project")

# 2. Increase download timeout & load packages
options(timeout = 600)
library(maftools)
library(TCGAbiolinks)
library(tidyverse)

# 3. Query and download TCGA-LIHC Somatic Mutations
query <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "MuTect2 Variant Aggregation and Masking"
)

GDCdownload(query)
lihc_raw <- GDCprepare(query)

# 4. Convert to a maftools MAF object
lihc_maf <- read.maf(maf = lihc_raw)

# 5. Load your candidate gene panel
gene_panel <- read_csv("data/gene_panel.csv")
target_genes <- gene_panel$`Gene Symbol`

# 6. Save Oncoplot figure to results/figures/
png("results/figures/circadian_oncoplot.png", width = 10, height = 8, units = "in", res = 300)
oncoplot(
  maf = lihc_maf,
  genes = target_genes,
  removeNonAltered = FALSE,
  draw_titv = TRUE,
  titleText = "Somatic Alterations in Circadian Gene Panel (TCGA-LIHC)"
)
dev.off()

# 7. Display Oncoplot directly in RStudio
oncoplot(
  maf = lihc_maf,
  genes = target_genes,
  removeNonAltered = FALSE,
  draw_titv = TRUE,
  titleText = "Somatic Alterations in Circadian Gene Panel (TCGA-LIHC)"
)