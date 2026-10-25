#======================================================================
# Author:      Farid
# Date:        July 2026
# Purpose: Merging files content in chronological order one after another
# Details: The following code are required to execute after extracting gene information with SNP information once gene mining is done based on SNP position. This will just combine the information chronologically, so that gene can be tracked from which SNP it was found.
#======================================================================
## Necessary libraries
library(tidyverse)
library(writexl)

# Step 1: Define the chronological structure of the files
# This matches the exact layout of groups per chromosome
chr_layout <- list(
  "1" = 1:6,   # Chr1 has group1 to group6
  "2" = 1:5,   # Chr2 has group1 to group5
  "3" = 1:10,  # Chr3 has group1 to group10
  "4" = 1:6,   # Chr4 has group1 to group6
  "5" = 1:12   # Chr5 has group1 to group12
)

# Step 2: Programmatically generate the filenames in strict chronological order
ordered_filenames <- c()

for (chr in names(chr_layout)) {
  for (grp in chr_layout[[chr]]) {
    # Construct the exact filename pattern
    file_name <- paste0("snps_chr", chr, "_group", grp, "_genelist_extracted_info.txt")
    ordered_filenames <- c(ordered_filenames, file_name)
  }
}

# Step 3: Verify all 39 files exist in the current folder before proceeding
missing_files <- ordered_filenames[!file.exists(ordered_filenames)]

if (length(missing_files) > 0) {
  warning("The following files were not found in the current directory:\n", 
          paste(missing_files, collapse = "\n"))
}

# Only proceed with files that actually exist
files_to_read <- ordered_filenames[file.exists(ordered_filenames)]

# Step 4: Import and merge the file contents sequentially
# Since your script outputted raw text without column headers, 
# we manually assign tidy headers here.
list_dataframes <- lapply(files_to_read, function(x) {
  read_tsv(x, 
           col_names = c("Pos", "Chr", "Gene_Start", "Gene_End", "Gene_ID"),
           col_types = "cciic") # c = character, i = integer
})
  
# Step 5: Stack all the individual dataframes into one master table
combined_dataframe<- dplyr::bind_rows(list_dataframes)

# Step 6: Save the perfectly ordered master data frame as an Excel spreadsheet
write_xlsx(combined_dataframe, path = "final_combined_geneid_list_38peaks.xlsx")

cat("Success! Combined", length(files_to_read), "files chronologically into 'merged_snps_master_list.xlsx'.\n")


#===================================================================
#Merging information from another file to this newly created combined #file using common cols
#===================================================================

## necessary packages
library(dplyr)
library(readxl)  
library(writexl)

## load master file, where I will add the necessary inf. from the another reference file
master_df <- read_xlsx("final_combined_geneid_list_38peaks.xlsx")

## Load the reference file containing the SNP and Group_No info
snp_reference <- read_xlsx("GWAS_SNP_information.xlsx")

## Clean and prepare the reference data
snp_ref_clean <- snp_reference %>%
  # Convert the raw numbers (1-5) to strings (Chr1-Chr5) right in the Chr column
  mutate(Chr = paste0("Chr", Chr)) %>%
  # Select only the keys and the new columns you want to bring over
  # (Make sure "Group_No" matches the exact spelling/symbols in your Excel header)
  select(Chr, Pos, SNP, Group_No)

## Perform the left_join for merging inf. from the reference file
final_master_data <- master_df %>%
  # Ensure the master Pos column is a number so it matches the reference file
  mutate(Pos = as.numeric(Pos)) %>% 
  # Join using the exact matching column names from both files
  left_join(snp_ref_clean, by = c("Chr", "Pos"))

## Save the final dataset to Excel
write_xlsx(final_master_data, "final_combined_geneid_list_38peaks1.xlsx")
