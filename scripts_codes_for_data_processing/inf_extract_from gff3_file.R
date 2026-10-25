library(tximport)
library(dplyr)
library(stringr)
library(writexl)
library(rtracklayer)

##importing
path_f<- file.path(".","Zmays_493_RefGen_V4.gene_exons.gff3")## giving file path
f_gff3<-rtracklayer::import(path_f)
head(f_gff3)

## converting gff3 file to dataframe to extract the target inf.
f_gff3_df<- as.data.frame(f_gff3)
head(f_gff3_df)

## extracting chr, type (gene) and gene name from gff3 file
inf_extract<- f_gff3_df %>%
  filter(type=="gene") %>%
  select(seqnames, Name)
  
## saving the result in excel file
write_xlsx(inf_extract, "gene_with_chr_list_maize_v4_annot.xlsx")
