#=======================================================================
# Author:      Farid
# Date:        July 2026
# Purpose: Subsetting and collapsing rows for each SNP.
# Details: Based on SNP, SNPs of peak GWAS region were subset and then, to make organize for manuscript,rows with different inf. in different columns were collapsed for same SNP. It can be used in such situation, where I need a single inf. in a single row with multiple inf. in a column by adding the inf with ";" or any other suitable separators
#=======================================================================

## necessary packages
library(readxl)
library(dplyr)
library(writexl)

# 1. Load your Excel files
master_data <- read_excel("Anox_combined_final_geneID_annt.xlsx")
snp_list    <- read_excel("anox_Chr3_top100_SNP.xlsx", sheet ="Sheet4")

# 2. Prepare the snp_list file 
# Keep only 'SNP' and 'Group_No.' (drops Model_name, Chromosome, and Position)
snp_list_clean <- snp_list %>%
  select(SNP, `Group_No.`)

# 3. Process, left join, and collapse the datasets
final_result <- snp_list_clean %>%
  # Left join preserves ALL rows from snp_list_clean
  left_join(
    master_data %>% select(-GeneID, -Annotation, -Group), 
    by = "SNP"
  ) %>%
  
  # Group by all fixed traits so they stay intact as individual columns
  group_by(SNP, `Group_No.`, Chr, Pos, Ref, Alt, MAF) %>%
  
  # Collapse the multi-row model statistics using semicolons
  # na.omit() ensures missing matches show up as clean, blank cells in Excel
  summarise(
    Model         = paste(na.omit(Model), collapse = "; "),
    P.value       = paste(na.omit(P.value), collapse = "; "),
    `H&B.P.Value` = paste(na.omit(`H&B.P.Value`), collapse = "; "),
    Effect        = paste(na.omit(Effect), collapse = "; "),
    .groups       = "drop"
  )

# 4. Preview the resulting dataset
print(final_result)

# 5. Export the final table to a new Excel file
write_xlsx(final_result, "anox_chr3_top100SNP_collapsed_output.xlsx")
