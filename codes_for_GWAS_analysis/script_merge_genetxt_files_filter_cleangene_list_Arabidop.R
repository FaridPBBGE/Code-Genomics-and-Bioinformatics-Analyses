# ===========================================================================
# Script Name: script_merge_genetxt_files_filter_cleangene_list
# Author:      Farid
# Date:        May 2026
# Purpose:     After using 4 bash scripts for extracting gene list based on GWAS SNP position, txt file for each peak region with having gene list need to be combined for downstream analysis. This script will provide final gene list for GWAS SNP. 

# ============================================================================

## Load required libraries
library(dplyr)
library(readr)
library(stringr)
library(writexl)

## List all files
files <- list.files(".", 
                    pattern = "\\.txt$", 
                    full.names = TRUE)

## Function to read and process each file
gene_merging_file <- function(file) {
  # remove suffix words from the file name
  file_tag <- str_remove(basename(file), "_gene_description\\.txt")

  # read file and add col names as no header in the original file
  df <- read.delim(file, 
                   header = FALSE, 
                   sep = "\t", 
                   col.names = c("transcript", "type", "description"), 
                   fill = TRUE, 
                   stringsAsFactors = FALSE)
  # add col. to track back data sources file
  df <- df %>% mutate(source = file_tag)
  return(df)
}

## Merge all files following newly created function ("gene_merging_file")
merged <- files %>% 
  lapply(gene_merging_file) %>% 
  bind_rows()

## Extract chromosome and group numbers for sorting to arrange the data chronologically
merged <- merged %>%
  mutate(
    chr = as.integer(str_extract(source, "(?<=chr)\\d+")),
    group = as.integer(str_extract(source, "(?<=group)\\d+"))
  ) %>%
  arrange(chr, group)

## remove decimal value from the gene IDs (because those are transcripts)
merged <- merged %>%
  mutate(gene = str_remove(gene, "\\.\\d+$"))

## check duplicates and counts
dup_summary<- merged %>%
  group_by(transcript) %>%
  summarise(count = n(), sources = paste(unique(source), collapse = ", "), .groups = "drop") %>%
  filter(count > 1)

## Remove duplicates (keep unique genes)
merged_unique <- merged %>%
  distinct(gene, .keep_all = TRUE)

## subset genes with unknown/blank descriptions and clean genes
unknown_transcript <- merged %>%
  mutate(description = str_trim(description)) %>%   # remove leading/trailing spaces
  filter(is.na(description) | description == "" | 
           str_detect(description, regex("\\bunknown function\\b", ignore_case = TRUE)))

clean_transcript <- merged %>%
  mutate(description = str_trim(description)) %>%   # same cleaning
  filter(!(is.na(description) | description == "" | 
             str_detect(description, regex("\\bunknown function\\b", ignore_case = TRUE))))

## save the results in excel format
write_xlsx(clean_transcript, "merged_clean_transcripts_gwas_peaks.xlsx")
write_xlsx(unknown_transcript, "excluded_unknown_transcripts_gwas_peaks.xlsx")

