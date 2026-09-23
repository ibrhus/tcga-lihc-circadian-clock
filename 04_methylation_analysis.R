# ==============================================================================
# Script: scripts/04_methylation_analysis.R
# Purpose: TCGA-LIHC Promoter DNA Methylation vs. mRNA Expression Correlation
# Project: Integrated Genomic & Clinical Analysis of Circadian Clock in HCC
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
cat("--> Loading required packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c("TCGAbiolinks", "SummarizedExperiment",
                   "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                   "ggplot2", "dplyr", "ggpubr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# 2. Query & Download TCGA-LIHC DNA Methylation (450k) Data
# ------------------------------------------------------------------------------
cat("--> Querying GDC for TCGA-LIHC DNA Methylation 450k data...\n")

query_meth <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value"
)

# Download files to 'GDCdata/' folder
cat("--> Downloading methylation files (this may take a few minutes)... \n")
GDCdownload(query_meth, method = "api", files.per.chunk = 10)

# Prepare SummarizedExperiment object
cat("--> Preparing SummarizedExperiment object...\n")
meth_se <- GDCprepare(query_meth)

# 3. Identify Promoter Probes for Target Genes
# ------------------------------------------------------------------------------
cat("--> Annotating promoter probes for target genes...\n")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

target_genes <- c("PER1", "CRY2", "NR1I3", "RORC")

# Filter probes in TSS1500, TSS200, 5'UTR, or 1stExon
promoter_probes <- as.data.frame(anno) %>%
  filter(UCSC_RefGene_Name %in% target_genes) %>%
  filter(grepl("TSS1500|TSS200|5'UTR|1stExon", UCSC_RefGene_Group)) %>%
  select(Name, chr, pos, UCSC_RefGene_Name, UCSC_RefGene_Group, Islands_Name, Relation_to_Island)

cat(sprintf("--> Found %d promoter probes across target genes.\n", nrow(promoter_probes)))

# 4. Extract Beta Values & Save Processed Checkpoint
# ------------------------------------------------------------------------------
beta_matrix <- assay(meth_se)[rownames(meth_se) %in% promoter_probes$Name, ]

# Create directory if it doesn't exist
if(!dir.exists("processed")) dir.create("processed")

saveRDS(list(beta = beta_matrix, probes = promoter_probes), "processed/meth_promoter_data.rds")
cat("--> Done! Checkpoint saved to 'processed/meth_promoter_data.rds'\n")