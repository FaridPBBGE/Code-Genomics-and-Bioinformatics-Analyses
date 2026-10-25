## Loading necessary packages for this scripts
library(writexl)
library(dplyr)
library(readr)

################## Processing GAPIT files #######################

### Making files with significant snps list and without maf0.05

## Defining the dir containing files
dir<- getwd()

## Getting a list of csv files in the dir
file_list<- list.files(path = dir, pattern = "*.csv", full.names = TRUE)

## Looping through each file
for (f in file_list) { 
  # read data from the file
  df<- read_csv(f)
  
  # Find the first p-value where BH value <=0.05
  threshold_data<- df[df$`H&B.P.Value` <= 0.05, ]
  if(nrow(threshold_data) > 0) {
    first_p_value<- threshold_data$P.value[1]
    # # Print the first P.value to the console
    print(first_p_value)
    # Print the message to the console
    cat(paste(f, "threshold level\n"))
  } else {
    cat(paste(f, "NO P.value meets the threshold level\n"))
  }
  # subset data where H.B.P.Value <=0.05
  subset_BHPvalue<- df[df$`H&B.P.Value` <= 0.05, ]
  # Check if subset_hbp is not empty before saving
  if (nrow(subset_BHPvalue) > 0) {
    # Save the subset data to a new Excel file
    write_xlsx(subset_BHPvalue, paste0("subset_5%BHP_", basename(f), ".xlsx"))
    # Subset the previously filtered data where MAF > 0.05
    subset_maf<- subset_BHPvalue[subset_BHPvalue$MAF > 0.05, ]
    
    # Check if subset_maf is not empty before saving
    if (nrow(subset_maf) > 0) {
      # Save the final subset data to another Excel file
      write_xlsx(subset_maf, paste0("subset_maf_", basename(f), ".xlsx"))
    } else {
      cat(paste(f, "No data meets the MAF  > 0.05 condition\n"))
    }
  } else {
  cat(paste(f, "No data meets the BHP value <=0.05 condition\n"))
  }
}


## Modifying the file for manhattan plot visualization in CMplot
library(dplyr)
library(readr)
dir_man<- getwd()

file_list<- list.files(path = dir_man, pattern = "*.csv", full.names = TRUE)
num_files <- length(file_list)
print(paste("Number of files loaded:", num_files))


for (file in file_list) {
  
  df_f_man<- read_csv(file) %>%
    select(chr,rs,ps,p_lrt) %>%
    rename(
      SNP=rs,
      Chromosome=chr,
      Position=ps,
      GEMMA=p_lrt
    ) %>%
    select(SNP, Chromosome, Position, GEMMA)

# Saving the renamed file
write_csv(df_f_man, paste0("Manhat_QQ_", basename(file), ".csv"))
}


###################### Processing GEMMA files ##########################

##### Calculating p-adjusted value using BH method from stats
library(dplyr)
dir<- getwd()

## Listing the files
file_list<- list.files(path = dir, pattern = "\\.txt$", full.names = TRUE)
num_files <- length(file_list)
print(paste("Number of files loaded:", num_files))


## Selecting necessary columns from the GEMMA result files
for (file in file_list) {
  # reading the file
  data_selected_col<- read.table(file, header = TRUE, sep = "\t") %>%
  # Selecting certain cols
  select(chr, rs, ps,allele1,allele0,af, beta, se, p_lrt)
  # Saving the output files
  write.table(data_selected_col, paste0("sel_cols_", basename(file)), sep="\t", row.names=FALSE, quote = FALSE)
}


## Calculating adjusted p-value for each file

for (file in file_list) {
  # reading the file
  data<- read.table(file, header = TRUE, sep = "\t")
  # Calculating adjusted pvalue
  data$adj_pvalue<- p.adjust(data$p_lrt, method = "BH")
  #saving the output file in csv format
  write.csv(data, paste0("Added_adjpvalue_", basename(file), ".csv"))
}


#### Making files with significant snps list and without maf0.05

dir<- getwd()

## Getting a list of csv files in the dir
file_list<- list.files(path = dir, pattern = "\\.txt$", full.names = TRUE)

## Looping through each file
for (file in file_list) {
  # read data from the file
  df<- read.table(file, header = TRUE, sep = "\t")
  
  # Find the first p-value where BH value <=0.05
  threshold_data<- df[df$`H&B.P.Value` <= 0.05, ]
  if(nrow(threshold_data) > 0) {
    first_p_value<- threshold_data$P.value[1]
    # # Print the first P.value to the console
    print(first_p_value)
    # Print the message to the console
    cat(paste(f, "threshold level\n"))
  } else {
    cat(paste(f, "NO P.value meets the threshold level\n"))
  }
  # subset data where H.B.P.Value <=0.05
  subset_BHPvalue<- df[df$`H&B.P.Value` <= 0.05, ]
  # Check if subset_hbp is not empty before saving
  if (nrow(subset_BHPvalue) > 0) {
    # Save the subset data to a new Excel file
    write_xlsx(subset_BHPvalue, paste0("subset_5%BHP_", basename(f), ".xlsx"))
    # Subset the previously filtered data where MAF > 0.05
    subset_maf<- subset_BHPvalue[subset_BHPvalue$MAF > 0.05, ]
    
    # Check if subset_maf is not empty before saving
    if (nrow(subset_maf) > 0) {
      # Save the final subset data to another Excel file
      write_xlsx(subset_maf, paste0("subset_maf_", basename(f), ".xlsx"))
    } else {
      cat(paste(f, "No data meets the MAF  > 0.05 condition\n"))
    }
  } else {
    cat(paste(f, "No data meets the BHP value <=0.05 condition\n"))
  }
}



## Loop through each file for calculating adjusted p-value
for (file in file_list) {
  # reading the file
  data<- read.table(file, header = TRUE, sep = "\t")
  # Calculating adjusted p-value using BH method
  data$p_adjusted<- p.adjust(data$p_lrt, method = "BH")
  # Naming output file name
  output_file<- paste0("Added_padj_", basename(file), ".csv")
  # Saving the output file in csv format
  write.csv(data, output_file)
}




######### Combining multiple files generated by different tools based on SNP #################

## loading packages
library(dplyr)
library(readxl)
library(readr)


## reading files
df_cmlm<- read_csv("Manhat_QQ_GAPIT.Association.GWAS_Results.CMLM.score.csv.csv") %>%
  select(SNP, Chromosome, Position, CMLM)

df_bl<- read_csv("Manhat_QQ_GAPIT.Association.GWAS_Results.BLINK.score.csv.csv") %>%
  select(SNP,Chromosome, Position, BLINK)

## reading MLMM files
df_ml<- read_csv("Manhat_QQ_GAPIT.Association.GWAS_Results.MLMM.score.csv.csv") %>%
  select(SNP, Chromosome, Position, MLMM)

## reading gemma file
df_gm<- read_csv("Manhat_QQ_Added_adjpvalue_sel_cols_AD_Anox.assoc.txt.csv.csv") %>%
  select(SNP, Chromosome, Position, GEMMA)


## performing full_join to combine all rows from all the files being combined
combined_file<- df_cmlm %>%
  full_join(df_bl, by = c("SNP","Chromosome", "Position")) %>%
  full_join(df_ml, by= c("SNP","Chromosome", "Position")) %>%
  full_join(df_gm, by = c("SNP","Chromosome", "Position")) 

##  Checking the columns order in the desired order
combined_file<- combined_file %>%
  select(SNP, Chromosome, Position, CMLM, BLINK,GEMMA,MLMM)

## Saving the final combined file
write_csv(combined_file, "AD__allsnpt_Anox_manhat_QQ_view.csv")

#### Merging files when each file has different column names  ###############
library(dplyr)
library(readxl)
library(writexl)

## loading data
f1<- read_excel("Anox_top1000snp_BLINK.xlsx")
f2<- read_excel("Anox_top1000snp_CMLM.xlsx")
f3<- read_excel("Anox_top1000snp_MLMM.xlsx")

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
write_xlsx(final_merged_file, "Anox_merged_all_models_top1000SNPs.xlsx")

############# Subsetting SNP data from the big result file ##################

## Loading data into R 
library(writexl)
library(readr)
library(dplyr)
dir<- getwd()

## Listing files
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

### Creating top 100/1000 SNPs based on top p-value of model and highlight them in the plot + For gene annotation analysis
library(dplyr)
library(writexl)
library(readxl)

data<- read_excel("AD_Mg_manhat_QQ_view.xlsx")

top_cmlm<- data %>% 
  arrange(CMLM) %>%
  slice(1:1000) %>%
  select(SNP,Chromosome, Position)

top_bl<- data %>% 
  arrange(BLINK) %>%
  slice(1:1000) %>%
  select(SNP,Chromosome, Position)

##top_gm<- data %>% 
  arrange(GEMMA) %>%
  slice(1:100) %>%
  select(SNP,Chromosome, Position)

top_ml<- data %>% 
  arrange(MLMM) %>%
  slice(1:1000) %>%
  select(SNP,Chromosome, Position)

## Saving the result
write_xlsx(top_cmlm, "Mg_top1000snp_CMLM.xlsx")
write_xlsx(top_bl, "Mg_top1000snp_BLINK.xlsx")
write_xlsx(top_ml, "Mg_top1000snp_MLMM.xlsx")


## saving the result with different column name for another software/analysis 
df_snp_cmlm<- data.frame(CMLM=top_cmlm)
df_snp_bl<- data.frame(BLINK=top_bl)
df_snp_ml<- data.frame(MLMM=top_ml)
##df_snp_gm<- data.frame(GEMMA=top_gm)

write_xlsx(df_snp_cmlm, "Cu_top1000snp_CMLM.xlsx")
write_xlsx(df_snp_bl, "Cu_top1000snp_BLINK.xlsx")
write_xlsx(df_snp_ml, "Cu_top1000snp_MLMM.xlsx")
##write_xlsx(df_snp_gm, "Anox_top100snp_GEMMA.xlsx")





####### Summarizing all genes information based on models to create ven diagram or summarize all the information

library(readxl)
library(dplyr)
library(openxlsx)
library(writexl)
library(tidyr)

##### Extracting distinct genes from all the genes found by various models (returns every distinct item, regardless of how many times they appear)

df<- read_excel("Anox_all_genes_by_all_models.xlsx")

## Extracting distinct genes
dist_genes<- df %>%
  distinct(Genes, .keep_all = TRUE)

## Saving the result
write_xlsx(dist_genes,"Anox_all_distinct_genes_by_all_models.xlsx")



########### Making a list of genes identify by how many models
## loading data
df<- read_excel("Anox_all_genes_by_all_models.xlsx")

## Find genes in all models/traits
genes_all_models<- df %>%
  group_by(Genes) %>%
  summarise(
    Models = n_distinct(Model),
    Models_names= paste(unique(Model),collapse = ",")) %>%
  filter(Models == length(unique(df$Model))) %>%
  select(Genes, Models_names)

## genes found in more than one model but not all models
genes_multiple_models<- df %>%
  group_by(Genes) %>%
  summarise(
    Models = n_distinct(Model),
    Models_names= paste(unique(Model),collapse = ",")) %>%
  filter(Models > 1 & Models < length(unique(df$Model))) %>%
  select(Genes, Models_names)

## Genes in only one model
genes_one_model<- df %>%
  group_by(Genes) %>%
  summarise(
    Models = n_distinct(Model),
    Models_names= paste(unique(Model),collapse = ",")) %>%
  filter(Models == 1) %>%
  select(Genes, Models_names)

## Combining all the result together
all_genes<- bind_rows(
  genes_all_models %>% mutate(Catagory= "All_Models"),
  genes_multiple_models %>% mutate(Catagory= "Multiple_Models"),
  genes_one_model %>% mutate(Catagory= "One_Model")
) %>%
  arrange(Genes)

## counting the catagory

all_genes<- all_genes %>%
  group_by(Catagory) %>%
  mutate(Catagory_count=  n()) %>%
  ungroup()



## making a list of gene under each catagory
list_gene_modelwise<- all_genes %>%
  group_by(Catagory) %>%
  summarise(genes_list=list(Genes)) %>%
  unnest(cols=c(genes_list)) %>%
  ungroup()

## Saving two different result object in two different sheet in the same excel file
wb<- createWorkbook()

# Add a sheet for the first result object
addWorksheet(wb, "gene-based_model_list")
writeData(wb, "gene-based_model_list", all_genes)

# for the second result object
addWorksheet(wb, "model-based_gene_list")
writeData(wb, "model-based_gene_list", list_gene_modelwise)

# saving both sheet together
saveWorkbook(wb, "Anox_summary_counting.xlsx", overwrite = TRUE)


## Creating venn diagram
library(VennDiagram)
library(readxl)
library(dplyr)

## Creating list as input file for venndiagram package
df<- read_excel("Anox_all_genes_by_all_models.xlsx")

genes_list<- list(
  CMLM= unique(df$Genes[df$Model =="CMLM"]),
  BLINK=unique(df$Genes[df$Model =="BLINK"]),
  MLMM=unique(df$Genes[df$Model == "MLMM"])
)

## Making venn diagram
ven_plot<- venn.diagram(
  x=genes_list,
  catagory.names= names(genes_list),
  filename = NULL,
  output= TRUE,
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex= 1.5,
  cat.cex = 1.5,
  cat.col = c("red", "blue", "green"),
  cat.pos=0
)

## saving the venn plot
png("Anox_venn_diagram.png")
grid.draw(ven_plot)
dev.off()

## Loading packages
library(dplyr)
library(readxl)

## Loading data
data<- read_excel("Zn_top1000snp_BLINK.xlsx")

trait<-"Zn"
model<- "BLINK"

for(chr in unique(data$BLINK.Chromosome)) {
  chr_data<- data %>% 
    filter(BLINK.Chromosome == chr) %>%
    select(BLINK.Position) %>%
    arrange(BLINK.Position)
  
  file_name<- paste0(trait,"_", model, "_chr", chr,"_top1000snp.txt")
  
  write.table(chr_data, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


#### Deleting GEMMA column from 50k subset SNP and creating file for each chromosome for chromosome wise manhattan plot visualization
with using function() to run multiple files
# Loading packages
library(dplyr)
library(readxl)
library(writexl)
library(tools)

### Function to process a single excel file
process_file<- function(file) {
  
## Getting base file name without extension
base_name_file<-file_path_sans_ext(basename(file))

## Loading data
df<- read_excel(file)

## Deleting GEMMA column
modify_df<- df %>%
  select(-GEMMA)

## saving the modified file
output_file_name<- paste0(base_name_file, "_noGEMMA.xlsx")
write_xlsx(modify_df, output_file_name)

## Subsetting newly modified data based on chromosome
chr_num<- unique(modify_df$Chromosome) %>% sort()

## Looping through each chromosome to subset and saving
for (chr in chr_num) {
  subset_modify_df<- modify_df %>%
    filter(Chromosome== chr)

# Defining file name for each subset to save
  file_name_subset<- paste0(base_name_file,"_Chr_", chr, ".xlsx")
  
# saving each subset
  write_xlsx(subset_modify_df, file_name_subset)
}
invisible(NULL)
}

## Using process_file function to multipe files

# List files, where function will be applied
list_files<- list.files(path = ".", pattern = "\\.xlsx$", full.names = TRUE)

# Applying the function to each file of the list
lapply(list_files, process_file)


##### Subsetting 1000SNPs from the GAPIT result file to obtain necessary cols information using looping

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


########### Merging two files based on two common cols ( case: merging geneId and annotation with snp file-GWAS-Arabidopsis) and order the cols as needed in the merged file using looping  and key term of the file name so that all files can be ran together

# Loading packages
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(writexl)


## matching between master and individual files

# reading files
  master_data <- read_csv("Zn_top1000SNP_snpeff_1K.ann_gwas_geneID_combined_id_ann.csv", show_col_types = FALSE)
  individual_data <- read_excel("Zn_1000SNP_inf_MLMM.xlsx")

# merging
  merged_data <- left_join(individual_data, master_data, by = c("Chr", "Pos"))
  
  final_col_order <- c(
    "Model", "SNP", "Chr", "Pos", "Ref", "Alt", "MAF", 
    "P.value", "H&B.P.Value", "Effect", "GeneID", "Annotation"
  )
  
  final_data <- merged_data %>%
    select(all_of(final_col_order))
  
  result_filename<- paste0("Zn_MLMM_final_geneID_annt.csv")
  write_xlsx(final_data, result_filename)
  
  
  
  
  
  





########### Adding a column at the last position with model name to track cols during merging files altogether using function and lapply

# Loading libraries
library(tidyverse)
library(readxl)
library(writexl)
library(tools)
## Loading files where new column will be added
list_files<- list.files(path = ".", pattern = "MLMM\\.xlsx$")

## Creating function for adding a column with filling up with a term for all the rows 
col_add<- function(file_name){
  
  # reading file
  df<- read_excel(file_name)
  
  # adding a column with filling up a certain term
  mod_df<- df %>%
    mutate(Model="BLINK")   ## If new col needs to relocate to a certain position, then use relocate(new colname,.after= use a interger to show after which-th col)
  
  #saving the result
  res_filename<- paste0(file_path_sans_ext(basename(file_name)),".xlsx")
  write_xlsx(mod_df, res_filename)
  
  # to see message for progress in console
  cat("Saved modified file:", res_filename, "\n")
  return(res_filename)
}

## applying the newly created function to all files using lapply
invisible(lapply(list_files, col_add))

