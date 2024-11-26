## Loading necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(writexl)

## Getting a list of filenames from the working directory that end in .tab/ any other particular pattern based on situation
f_names<- list.files(pattern ="*geneidReadsPerGene.out.tab") 

## Reading each file in a list of files, extracting only the first and second columns after skipping the first four rows
dfs<- lapply(f_names, function(file) {
  read_tsv(file, col_names = c("Gene", file),skip = 4, col_types = cols_only(Gene= col_character(), !!file := col_double()))
})

## Combining all dataframes into one
combined_df<-reduce(dfs, full_join, by="Gene")

## Removing duplicate gene columns
combined_df<- combined_df %>% distinct(Gene, .keep_all = TRUE)

## Saving the result in excel
write_xlsx(combined_df, "compiled_counts.xlsx")
