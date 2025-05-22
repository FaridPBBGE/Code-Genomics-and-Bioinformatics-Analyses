#####################Processing salmon_output_files ########################

## loading packages
library(tximport)
library(dplyr)
library(ggplot2)
library(stringr)
library(writexl)

## List all zip files in the current directory
f_zipped<- list.files(path = ".", pattern = "\\.zip$")

## loop through the list of files, make folder, name the folder and unzip file
for (zf in f_zipped) {
  # remove the .zip extension to create folder name
  folder_name<- gsub("\\.zip$","",zf)
  
  # create folder using folder name jsut created now
  dir.create(folder_name, showWarnings = FALSE)
  cat(paste("Unzipping", zf, "to folder", folder_name, "\n"))
  
  # unzip the file into newly created folder
  unzip(zipfile = zf, exdir = folder_name)
  cat(paste("Finished unzipping:", zf, "\n"))
}


## getting files salmon folders
c_files<- list.files(path =".", pattern = ".sf", full.names = TRUE, recursive = TRUE)
## renaming files with accession number
names(c_files)<- stringr::str_split(c_files,pattern = "/", simplify = TRUE)[,2] %>% ##[] will select which column I need to select for replace. The column will be output of str_split code
stringr::str_replace("_Rquant","") ## removing quant 


## loading gtf file of reference genome using rtracklayer package
BiocManager::install("rtracklayer")
library(rtracklayer)
##importing
path_f<- file.path(".","Zmays_493_RefGen_V4.gene_exons.gff3")## giving file path
f_gff3<-rtracklayer::import(path_f)
head(f_gff3)

## converting gff3 file to dataframe to extract the target inf.
f_gff3_df<- as.data.frame(f_gff3)
head(f_gff3_df)
## getting transcript and gene ID from this df file, so that it can be used for salmon count files (converting transcript level count to gene-level count purpose)
t2gene<- f_gff3_df%>%
  dplyr::filter(type== "mRNA")%>% ### for gtf file, "mRNA" will be "transcript"
  select(ID,Parent)%>%
  mutate(ID=str_replace(ID,pattern= "\\.RefGen_V4$",""),
         Parent=str_replace(Parent,pattern="\\.RefGen_V4$",""))

## importing t2gene file into tximport
t2gene_test<- tximport(c_files,type = "salmon",tx2gene = t2gene)

## Extract gene-level tpm data (## It is for extracting tpm and readcounts for other use outside of Deseq2)
tpm_gene<- t2gene_test$abundance
head(tpm_gene, n=5)

## Before converting to dataframe, rownames need to extracted from tpm_gene and need to be added to dataframe after convertiong to dataframe becasue tpm_gene is matrix.
gene_names<- row.names(tpm_gene)

## Converting to dataframe
tpm_gene_df<- as.data.frame(tpm_gene)

## Adding gene names to the dataframe
tpm_gene_df$Gene<- gene_names
tpm_gene_df <- cbind(Gene = gene_names, tpm_gene_df) # Add gene names as the first column

## Saving into excel file
write_xlsx(tpm_gene_df,"TPM_Hi-maize-A_2025_Salmon.xlsx")


##################### Processing TPMCalculator_output_files #############################

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
