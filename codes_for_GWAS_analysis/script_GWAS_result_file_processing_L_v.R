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

## Calculating p-adjusted value using BH method from stats
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

### Making files with significant snps list and without maf0.05

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




######### Combining multiple files generated by different software based on SNP #################

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

data<- read_excel("AD_Zn_manhat_QQ_view.xlsx")

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

## saving the result
df_snp_cmlm<- data.frame(CMLM=top_cmlm)
df_snp_bl<- data.frame(BLINK=top_bl)
df_snp_ml<- data.frame(MLMM=top_ml)
##df_snp_gm<- data.frame(GEMMA=top_gm)

write_xlsx(df_snp_cmlm, "Zn_top1000snp_CMLM.xlsx")
write_xlsx(df_snp_bl, "Zn_top1000snp_BLINK.xlsx")
write_xlsx(df_snp_ml, "Zn_top1000snp_MLMM.xlsx")
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
