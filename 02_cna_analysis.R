# --------------------------------------------------
# 02_cna_analysis.R - Copy Number Alterations Analysis
# --------------------------------------------------

library(curatedTCGAData)
library(maftools)
library(tidyverse)

# 1. Load gene panel
gene_panel <- read_csv("data/gene_panel.csv", show_col_types = FALSE)
target_genes <- gene_panel$`Gene Symbol`

# 2. Fetch GISTIC discrete thresholded CNA data for LIHC
lihc_cna <- curatedTCGAData(
  diseaseCode = "LIHC",
  assays = "GISTIC_ThresholdedByGene",
  version = "2.1.1",
  dry.run = FALSE
)

# 3. Extract matrix and set row names from Gene.Symbol metadata
se_cna <- lihc_cna[[1]]
cna_matrix <- assays(se_cna)[[1]]
rownames(cna_matrix) <- rowData(se_cna)$Gene.Symbol

# 4. Filter for circadian target genes
circadian_cna <- cna_matrix[rownames(cna_matrix) %in% target_genes, , drop = FALSE]

# 5. Format long table with matched 15-character TCGA barcodes
cna_long <- circadian_cna %>%
  as.data.frame() %>%
  rownames_to_column(var = "Hugo_Symbol") %>%
  pivot_longer(-Hugo_Symbol, names_to = "Tumor_Sample_Barcode", values_to = "CNA_Status") %>%
  filter(CNA_Status %in% c(-2, 2)) %>%
  mutate(
    CN = ifelse(CNA_Status == 2, "Amp", "Del"),
    Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 15)
  ) %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode, .keep_all = TRUE) %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, CN)

# 6. Save processed CNA table
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write_csv(cna_long, "data/processed/circadian_cna_maftools.csv")

# 7. Merge MAF + CNA data
lihc_maf_cna <- read.maf(
  maf = "data/raw_maf/data_mutations.txt",
  cnTable = cna_long
)

# 8. Save combined Oncoplot
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

png("results/figures/circadian_oncoplot_with_cna.png", width = 11, height = 8.5, units = "in", res = 300)
oncoplot(
  maf = lihc_maf_cna,
  genes = target_genes,
  draw_titv = TRUE,
  titleText = "Somatic Alterations & CNA in Circadian Panel (TCGA-LIHC)"
)
dev.off()