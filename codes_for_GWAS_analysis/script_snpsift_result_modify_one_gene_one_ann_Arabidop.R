# ============================================================================
# Script Name: script_snpsift_result_mod_one_gene_one_ann.R
# Author:      Farid
# Date:        May 2026
# Purpose:     This script will be used after applying bash scripts for two                                 steps modification to modify the original snpsift-derevied vcf                               file. As a third step, this R-script will be used to modify     snfsift                      -derived result table to make one gene-one annotation-one row with counting                  column to see for how many times that gene has been hit by SNPs. 
# ============================================================================


## Load required library
library(dplyr)
library(openxlsx)

summary_df <- read.delim("snpeff_38gwas_peaks_genes_ann_for_R_final.txt", sep="\t", header=TRUE, check.names=FALSE) %>%
  group_by(`ANN[0].GENEID`) %>%
  summarise(
    EFFECTS     = paste(sort(unique(`ANN[0].EFFECT`)), collapse = ", "),
    IMPACTS     = paste(sort(unique(`ANN[0].IMPACT`)), collapse = ", "),
    TOTAL_COUNT = n()
  ) %>%
  arrange(desc(TOTAL_COUNT))


write.xlsx(summary_df, "38gwas_peaks_snpeff_gene_ann_summary.xlsx")
