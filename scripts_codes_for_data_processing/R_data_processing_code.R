##------------------------------------------------------------------------------

#Converting excel file from csv to excel files for multiple files

##------------------------------------------------------------------------------
## Loading libraries
library(openxlsx)
library(readr)

input_folder <- "."  # change if needed

# List all csv files in the folder
csv_files <- list.files(path = input_folder, pattern = "\\.csv$", full.names = TRUE)

# Loop through each XLSX file and save as CSV
for (file in csv_files) {
  
  # Read csv
  data <- read_csv(file, col_names =TRUE)
  
  # Construct file name
  excel_file <- paste0(tools::file_path_sans_ext(file), ".xlsx")
  
  # Write as CSV
  write.xlsx(data, excel_file)
  
  cat("Converted:", basename(file), "->", basename(excel_file), "\n")
}










########### Purpose: Combine excel files which has common columns and there are 7 types of files in terms of name and each term has 3 more sub-types-Meging files using looping for automation

## lOading required packages
library(dplyr)
library(readxl)
library(writexl)
library(purrr)

## getting files list
all_files_list<- list.files(".", pattern = ".*_final_geneID_annt\\.xlsx$")

## finding key terms so that files can be combined based on sub-key terms
key_terms<- unique(sub("^(.+?)_.*$", "\\1", all_files_list))

## looping through each key term for combining 3 sub-key terms files
for (term in key_terms) {
  file_list<- all_files_list[grepl(paste0("^", term, "_"), all_files_list)]
  
  # reading, sorting, and combining 3 sub-key term files
  combined_files<-lapply(file_list, function(file) {
    
    #reading
    read_excel(file) %>%
      
      # ensuring numeric columns are treated as numbers
      mutate(
        `P.value` = as.numeric(`P.value`),
        `H&B.P.Value` = as.numeric(`H&B.P.Value`),
        `Effect` = as.numeric(`Effect`),
        `MAF` = as.numeric(`MAF`),
        # Also ensure Chr and Pos are treated as numeric for consistency
        Chr = as.numeric(Chr),
        Pos = as.numeric(Pos)
      ) %>%
      
      #sorting
      arrange(P.value)
  }) %>%
    bind_rows()
  
  # creating new file name for the combined file
  filename_output<- paste0(term, "_combined_final_geneID_annt.xlsx")
  
  # saving the combined file in excel format
  write_xlsx(combined_files, filename_output)
  
  cat("Combined files for:", term, "into", filename_output, "\n")
  }
cat("All files have been processed and combined")

### merging all combined files into a single excel file with different sheet using key term of the files

## getting files list
all_combined_files<- list.files(pattern =".*combined_final_geneID_annt\\.xlsx$")

## extract key terms from the file names
key_terms<- unique(sub("^(.+?)_.*$", "\\1", all_combined_files))

## creating names list to make sure ordering is correct
file_list_with_names<- setNames(all_combined_files, key_terms)

## reading each file to merge into a single excel file
list_dataframes<- map(file_list_with_names, read_excel)

## saving the combine files into a single excel file
write_xlsx(list_dataframes, "Final_GWAS_all_traits_top1000SNP_gene_annt.xlsx")

###### XXXXDoneXXXXX #################


############ Purpose: Finding to extract top 30 and bottom 30 samples for a list of variables using "function" to use it in loop  ###############

## loadign necessary library
library(readxl)
library(dplyr)

# Reading the data
df<- read_excel("Phenotype.xlsx")

# Variables need to be processed
variables<- c("Ca","Cu","Fe","Mg","Mn","Zn")

# Creating a function to extract top 30 and bottom 30 samples for a list of variables
extract_top_bottom<- function(data, var, n=30) {
  top<- data %>%
    arrange(desc(.data[[var]])) %>%
    slice_head(n=n) %>%
    select(Taxa)
  
  bottom<- data %>%
    arrange(.data[[var]]) %>%
    slice_head(n=n) %>%
    select(Taxa)
  
  result<- bind_rows(top, bottom) %>%
    mutate(Group= rep(c("Top","Bottom"), each=n),
           variable= var)
  
  return(result)
}

# Applying the function to my data
final_result<- lapply(variables, function(v) extract_top_bottom(df,v)) %>%
  bind_rows()

# Viewing the result
head(final_result)

# Savign the result
library(writexl)
write_xlsx(final_result, "top_bottom_30samples.xlsx")


######### Purpose: Extracting randomly 30 samples from categories data ##################

## loadign necessary library
library(readxl)
library(dplyr)

# Reading data
df_cat<- df[,c(1,8)]

# Extracting randomly 30 samples from two categories each
sampling_30<- df_cat %>%
  filter(score %in% c(1,3)) %>%
  group_by(score) %>%
  slice_sample(n=30) %>%
  ungroup() %>%
  select(Taxa, score)

# Savign the result
library(writexl)
write_xlsx(sampling_30, "anoxia_60samples.xlsx")

~~~~~~~~~~~~~~Purpose: Merging files~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   ### S-1: Adding rows from one data to another data where there is one common col between two data and other columns's name are different

## loadign necessary library
library(readxl)
library(dplyr)

# Reading data
data1<- read_excel("top_bottom_30samples.xlsx")
data2<- read_excel("anoxia_60samples.xlsx")

# modifying data1 to match columns between 2 data by adding empty column that is available in data2 only
data1_mod<- data1 %>%
  mutate(score=NA)

# modifying data2 to match columns between 2 data by adding empty columns that are available in data1 only
data2_mod<- data2 %>%
  mutate(Group=NA, variable=NA) %>%
  select(Taxa, Group, variable, score)

# binding them together
final_data<- bind_rows(data1_mod, data2_mod)

# Savign the result
library(writexl)
write_xlsx(final_data, "merged_60samples_all_traits.xlsx")

#### S-2: Merging files when each file has different column names
f1<- read_excel("Mg_top1000snp_BLINK.xlsx")
f2<- read_excel("Mg_top1000snp_CMLM.xlsx")
f3<- read_excel("Mg_top1000snp_MLMM.xlsx")

## Need to rename the cols of each file to make it same before merging them
f1<- f1 %>% rename(SNP= BLINK.SNP, Chr= BLINK.Chromosome, Pos= BLINK.Position)
f2<- f2 %>% rename(SNP=CMLM.SNP, Chr= CMLM.Chromosome, Pos= CMLM.Position)
f3<- f3 %>% rename(SNP= MLMM.SNP, Chr= MLMM.Chromosome, Pos= MLMM.Position)

## Merging files
merged_file<- f1 %>%
  full_join(f2, by= c("SNP","Chr", "Pos")) %>%
  full_join(f3,  by = c("SNP","Chr", "Pos"))

## Arranging the merged data based on Chr column
final_merged_file<- merged_file %>%
  arrange(Chr, Pos)

## Saving the merged file
write_xlsx(final_merged_file, "Mg_merged_all_models_top1000SNPs.xlsx")


####### Purpose: To extract unique value from a single column that has name    ################
library(readxl)
library(dplyr)

# reading data
df<- read_excel("merged_60samples_all_traits.xlsx")

# Method-1: using base R
unique_samples<- unique(unlist(df$Taxa))

# converting dataframe
unique_sample_df<- data.frame(Taxa_unique=unique_samples)

# Savign the result
library(writexl)
write_xlsx(unique_sample_df, "unique_taxa.xlsx")

# Method-2: using dplyr
uniq_sample_dplyr<- df %>%
  select(Taxa) %>%
  distinct()

# converting dataframe
uniq_sample_dplyr_df<- data.frame(Taxa_unique= uniq_sample_dplyr$Taxa)

# Savign the result
write_xlsx(uniq_sample_dplyr_df, "unique_taxa_dplyr.xlsx")


~~~~~~~~~~Purpose: Subsetting data from the master file~~~~~~~~~~~~~~
  
####S-1:

# using base R
library(readxl)
library(writexl)

# reading subset file
subset_list<- read_excel("mapped_marker_list.xlsx")            
# This is the list of names of variable that needs to be subset from the master file

# Reading master file
master_file<- read_excel("format_2_F2_double_letter_format.xlsx")

# Subsetting 
subset_data<- master_file[master_file[[1]] %in% subset_list[[1]], ]

# using dplyr
library(dplyr)
subset_data<- master_file %>%
  filter(.data[[1]] %in% subset_list[[1]])

# Subsetting specific number of columns
subset_data<- master_file %>%
  select(1:11) %>%
  filter(`rs#` %in% subset_list[[1]])

# saving result in excel
write_xlsx(subset_data, "1009_SNP_for_vcf.xlsx")


##### S-2: Subsetting 10,000 SNPs based on lowest p-value per chr for multiple files using loop
file_list<- list.files(path = dir, pattern = ".csv", full.names = TRUE)
num_files <- length(file_list)
print(paste("Number of files loaded:", num_files))

## Subsetting 10,000 SNPs based on lowest p-value per chr.
for (file in file_list) {
  # reading each file
  subset_50Ksnps<- read_csv(file) %>%
    group_by(Chromosome) %>%
    arrange(GEMMA) %>%
    slice_head(n=10000) %>%
# Sorting the subset data based on chr and snp's position
  arrange(Chromosome, Position)
  
# Saving the output as excel with new name
write_xlsx(subset_50Ksnps, paste0("subset50k_snp_", basename(file), ".xlsx"))
}

#### S-3: Selecting only 4 cols needed from a file 
## Modifying the subset file for manhattan visualization
df_file<- read_excel("process_GAPIT.Association.GWAS_Results.CMLM.Fe.xlsx")

## Selecting only 4 cols needed for CMplot for manhattan and QQ plot visualization
df_file_subset<- df_file %>%
  select(SNP, Chr,Pos,P.value) %>%
  rename(
    SNP=SNP,
    Chromosome=Chr,
    Position=Pos,
    CMLM=P.value
  )
head(df_file_subset)

## Saving the renamed file
write_xlsx(df_file_subset, "AD_Fe_for_man_QQ_plot.xlsx")


##### S-3: Subsetting using loop, where multiple master files and individual files are located and data need to be subset based on matched master and individual file
# Loading packages
library(tidyverse)
library(readxl)
library(writexl)
library(tools)

## Master file list (from which data will be subsetting)
master_files<- list.files(path = ".", pattern = "\\.csv$", full.names = TRUE)

##List of files with having list of snps (1000snps), which data will be used for subsetting
ind_files<- list.files(path = ".", pattern = "\\.xlsx$")

## Loop through each master file to match with corresponding individual file and subsetting data
for (mf in master_files) {
  # Finding key term in master file to match with individual files
  key_term<- str_extract(mf, "(?<=MLMM\\.)(.*?)(?=\\.csv$)")
  
  # Skip to the next file if no key term is found
  if (is.na(key_term)) {
    next
  }
  # Finding individual files based on key term found in master file
  matching_ind_file<- ind_files[str_detect(ind_files, key_term)] 
  
  # checking if a match was found before reading and subsetting tasks
  if (length(matching_ind_file) > 0) {
    cat("Processing file for:", key_term, "\n")
    
    # reading master and individual file
    mf_data<- read_csv(mf, show_col_types = FALSE)
    ind_files_data<- read_excel(matching_ind_file)
  
    # extract col from the individual file for which data will be subset from master file
    col_ind_file<- pull(ind_files_data, var= 1)
    
    # subsetting master file based extracted col from ind.file
    subset_data<- mf_data %>%
      filter(mf_data[[1]] %in% col_ind_file) %>%
      select(-6)
    
    # saving the result in excel file
    result_filname<- paste0( key_term, "_1000SNP_inf_MLMM.xlsx")
    result<- write_xlsx(subset_data, result_filname)
    
    cat("Successfully created:", result_filname, "\n\n")
  }
else {
  cat("No matching individual file found for:", mf, "\n")
}  
}


#### S-4: Based on a common column (SNP) between two datasets,the goal is to merge information from the source file (Model column in Data2) into the target file (Data1). Since multiple models can correspond to the same SNP, all matching models from Data2 should be collapsed and placed into a single cell in the merged Data1, ensuring each SNP in Data1 retains one row with its associated models listed together

## Loading required packages
library(readxl)
library(dplyr)
library(writexl)

## Reading the data
target_file<- read_excel("Anox_merged_all_models_top1000SNPs.xlsx", sheet = "Sheet1")
source_file<- read_excel("Anox_merged_all_models_top1000SNPs.xlsx", sheet = "Sheet2")

## Collapsing multiple models per SNP into one cell before merging with the target file
collapsed_source_file<- source_file %>%
  group_by(SNP) %>%
  summarise(Model=paste(unique(Model), collapse = ";"), .groups = "drop")

## Merging target col. inf. based on common col.
merged_file<- target_file %>%
  left_join(collapsed_source_file,by="SNP") %>%
  select (Model, SNP, Chr, Pos)

## Saving the result
write_xlsx(merged_file, "Anox_merged_all_models_top1000SNPs_1.xlsx")

  


## Compiled from another script file
library(readxl)
df<-read_excel("Book1.xlsx")

library(tidyverse)
## clean data
df_clean<- df %>%
  mutate(
    clean_sample_k=str_remove_all(sample_k,"\\W") %>% str_to_lower(),
    clean_sample_reg=str_remove_all(sample_reg,"\\W") %>% str_to_lower()
    )

## Create a logical vector for gettign matched names
df_clean<- df_clean %>%
  mutate(match=clean_sample_k %in% clean_sample_reg)

## Filtering common names
common_names<- df_clean %>%
  filter(match) %>%
  pull(clean_sample_k)

## Counting common names 
length(common_names)

## saving file
library(writexl)
df_common_names<- as.data.frame(common_names)
write_xlsx(df_common_names,"common_names.xlsx")
