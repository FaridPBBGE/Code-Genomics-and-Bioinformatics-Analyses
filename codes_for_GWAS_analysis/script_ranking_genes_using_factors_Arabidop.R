# =======================================================================
# Script Name: script_ranking_genes_using_factors.R
# Author: Farid
# Date:May 2026
# Purpose: Combine various cols with having various inforamtion (such gene expression, snpeff annotation and so on) based on geneid and try to rank genes
#Use: This script was used to rank genes for Arabidopsis gene for 38 GWAS peaks found in 5 Chr. based on top 100 SNP in each chromosome.

# ======================================================================

#---------------------------------------------------------------
#reshaping file2 (inf will be added to file1); main purpose is to find #which genes are found in which genotype; then merge this inf. to the #file1
#----------------------------------------------------------------

## load required packages
library(tidyverse)
library(readxl)
library(openxlsx)
library(writexl)

## read the file, where I will add new inf. from another file
final_df<- read_excel("scale_for_various_factors.xlsx", sheet = "unknown_gene_score_anox_regen")

## read new file, which inf. will be added to final_df file
new_inf_df<- read_excel("combined_DEG_wt_mt.xlsx")

## For each gene in file2, collect genotype presence
genotype_presence <- new_inf_df %>%
  group_by(gene) %>%
  summarise(
    has_WT   = any(genotype == "WT"),
    has_MT = any(genotype == "MT"),
    .groups = "drop"
  ) %>%
  mutate(
    genotype_presence = case_when(
      has_WT & has_MT ~ "Both",
      has_WT & !has_MT ~ "WT only",
      !has_WT & has_MT ~ "MT only",
      TRUE ~ "No match"
    )
  ) %>%
  select(gene, genotype_presence)

## Merge into the final file
merged_df <- final_df%>%
  left_join(genotype_presence, by = "gene")

## save the result in pre existing excel file
wb<- loadWorkbook("scale_for_various_factors.xlsx")

writeData(wb, sheet = "unknown_gene_score_anox_regen", x=merged_df)

saveWorkbook(wb, "scale_for_various_factors.xlsx", overwrite = TRUE)


#---------- creating new column based on a value of a col and then merging to the file1 with common col name --------

#### load required packages
library(dplyr)
library(readxl)
library(openxlsx)
library(writexl)

## read the file, where I will add new inf. from another file
final_df<- read_excel("merged_clean_transcripts_gwas_peaks.xlsx")

## read new file, which inf. will be added to final_df file
new_inf_df<- read_excel("gwas_clean_genes_ann.xlsx")

## keep only the common col + info col, and rename the info col
new_inf_df <- new_inf_df %>%
  select(gene, `HN-score`) %>%
  rename(hypoxia_meta_HN_score =`HN-score` )

## convert transcript IDs into single gene IDs (remove decimal part)
final_df <- final_df %>%
  mutate(gene = sub("\\..*", "", transcript)) # keep only part before the dot

## Collapse HN-scores rows per gene
gene_scores<- final_df %>%
  group_by(gene) %>%
  summarise(`HN-score` = paste(na.omit(`HN-score`), collapse = ","), .groups = "drop") %>%
  rename(hypoxia_meta_HN_score =`HN-score`)


## create new column "expression" based on sign of logFC_R_subm
new_inf_df<- new_inf_df %>%
mutate(shoot_exp_tolXtreat = if_else(logFC_S_subm > 0, "Up", "Down"))

## merge into main file
merged_df<- new_inf_df %>%
  left_join(gene_scores, by= "gene")

## save the result in pre existing excel file
wb<- loadWorkbook("merged_clean_transcripts_gwas_peaks.xlsx")

writeData(wb, sheet = "Sheet1", x=merged_df)

saveWorkbook(wb, "merged_clean_transcripts_gwas_peaks.xlsx", overwrite = TRUE)

# alternative way to update the file uisng openxl
write.xlsx(merged_df, "gwas_clean_genes_ann.xlsx", overwrite = TRUE)


#--------------------------------------------------------------------
#Merging GO ID and description into the meta file containing gene #list by #collapsing multiple rows for one gene into onerow for one gene #format
#--------------------------------------------------------------------

## load required packages
library(dplyr)
library(readxl)
library(openxlsx)
library(writexl)

## read GO annotation file
go_ann<- read_excel("clean_gene_go_ann.xlsx")

## filter for biological process information only from the file
bp_ann<- go_ann %>% filter(category == "proc")

## modifying for one row-one gene format (collapsing multiple rows for same gene's information)
collapsed_row<- bp_ann %>%
  group_by(Locus) %>%
  summarise(
    GO_ID = paste(unique(`GO ID`), collapse = ";"),
    GO_Term = paste(unique(`GO term`), collapse = ";")
  ) %>%
  ungroup()

## Save the result in excel
write_xlsx(collapsed_row, "collapsed_clean_gene_go_ann.xlsx")

## merge new inf. from this newly created go annotation file into the main file
#read the file, where I will add new inf. from another file
final_df<- read_excel("gwas_clean_genes_ann.xlsx")

## read new file, which inf. will be added to final_df file
new_inf_df<- read_excel("collapsed_clean_gene_go_ann.xlsx")

## merging to the final_df file
merged_inf<- final_df %>%
  left_join(new_inf_df, by = c("gene" = "Locus"))

## save the result
write.xlsx(merged_inf, "gwas_clean_genes_ann.xlsx", overwrite = TRUE)
