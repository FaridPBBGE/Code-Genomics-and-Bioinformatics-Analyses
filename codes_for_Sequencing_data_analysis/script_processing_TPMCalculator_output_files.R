library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(writexl)
library(tools)

## reading master file
master_gene_list<- read_excel("master_gene_list_v4.xlsx")

##
file_list<- list.files(pattern = "*_geneidAligned.sortedByCoord.out_genes_extracted.xlsx")

##
merged_dfs<- list()

## looping through each file and merge them with master file
for (f in file_list){
  # reading each excel file
  temp_df<-read_excel(f, col_names = TRUE)
  names(temp_df)<- c("Gene", "TPM")
  #doing merging by left_join to keep all rows from df1 and the matching rows from df2; non-matching rows from df2 will have NA
  merged_df<- left_join(master_gene_list, temp_df, by="Gene")
  #keeping only Gene and TPM cols from merged df
  merged_df<- select(merged_df, Gene, TPM)
  # Extrating file name without extension and path
  file_name<- tools::file_path_sans_ext(basename(f))
  #Renaming TPM col based on original filename
  colnames(merged_df)[2]<- paste0("TPM_",file_name)
  # adding merged df to list
  merged_dfs[[file_name]]<- merged_df
}

## combining all dfs into one by joining on the Gene col
combined_df<- reduce(merged_dfs, full_join, by="Gene")
print(head(combined_df))
## saving the result
write_xlsx(combined_df,"TPM_TPMCalculator_V4_maize_Hi-A.xlsx")
