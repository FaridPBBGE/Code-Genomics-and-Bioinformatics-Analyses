## loading DEVis package
library(DESeq2)

## Installing devtools for installing form github
if(!require("devtools")) install.packages("devtools")
library(devtools)

## installing DEvis
devtools::install_github("price0416/DEvis/DEVis")

## loading necessary packages
library(tximport)
library(dplyr)
library(ggplot2)
library(stringr)
library(writexl)

##### getting ready files for DEVis package
c_files<- list.files(path =".", pattern = ".sf", full.names = TRUE, recursive = TRUE)
## renaming files with accession number
names(c_files)<- stringr::str_split(c_files,pattern = "/", simplify = TRUE)[,2] %>% ##[] will select which column I need to select for replace. The column will be output of str_split code
  stringr::str_replace("_quant","") ## removing quant 


## loading gtf file of reference genome using rtracklayer package
BiocManager::install("rtracklayer")
library(rtracklayer)
##importing
path_f<- file.path(".","Zmays_RefGen_V4.gene_exons.gff3")## giving file path
f_gff3<-rtracklayer::import(path_f)
head(f_gff3)

## converting gff3 file to dataframe to extract the target inf.
f_gff3_df<- as.data.frame(f_gff3)
head(f_gff3_df)
## getting transcript and gene ID from this df file, so that it can be used for salmon count files (converting transcript level count to gene-level count purpose)
t2gene<- f_gff3_df%>%
  dplyr::filter(type== "mRNA")%>% ### for gtf file, mRNA will be transcript
  select(ID,Parent)%>%
  mutate(ID=str_replace(ID,pattern= "\\.RefGen_V4$",""),
         Parent=str_replace(Parent,pattern="\\.RefGen_V4$",""))

## importing t2gene file into tximport
t2gene_test<- tximport(c_files,type = "salmon",tx2gene = t2gene)

###### For multi-factor experiment############

### Creating target data (equivalent to meta file in DESeq2)
## Creating Exp. condition/factor 
target_data<-data.frame(condition=c("Early", "Early", "Early", "Early", "Late", "Late", "Late", "Late"), type=c("A","A","B","B","A","A","B","B"))

## naming the condition file's row name according to the columns of counts
rownames(target_data)<-colnames(t2gene_test$counts)

## Converting characters to factor
target_data$type<- factor(target_data$type)

## saving target data as csv file as required by DEVis package
write.csv(target_data, file = "target_matrix.csv", row.names =TRUE, col.names = TRUE)

### Creating count data for DEvis 
count_data_f<- as.matrix(t2gene_test$counts)
rownames(count_data_f)<- rownames(t2gene_test$counts)
colnames(count_data_f)<- colnames(t2gene_test$counts)
## Readjusting the cols names to make the first column headerless
count_data_f<-cbind(Genename=rownames(count_data_f),count_data_f)
row.names(count_data_f)<- NULL

## saving count data as txt file as required by DEVis package
write.table(count_data_f, file="count_matrix_1.txt", sep = "\t", row.names = FALSE, col.names = c("", colnames(t2gene_test$counts)), quote = FALSE)

                  #### Starting DEVis package #######
library(DEVis)
## Define the base dir for the analysis
base_dir<- getwd()

## Creating the dir structure for results
create_dir_struct(base_dir)

## Specifying input directories
cnt_dir<- "C:/Users/asmfaridul.islam/OneDrive - Texas A&M AgriLife/Documents/R/R_Output/CS_postdoc/Deseq2_practice/DEVis/counts/"

tgt_dir<- "C:/Users/asmfaridul.islam/OneDrive - Texas A&M AgriLife/Documents/R/R_Output/CS_postdoc/Deseq2_practice/DEVis/targets/"

init_data_paths(cnt_dir, tgt_dir)

## Specifying count and target file names
count_matrix<- "count_matrix_1.txt"
target_matrix<- "target_matrix.csv"

## Reading count and target data
count_data<- prep_counts(count_matrix)
## making interger to run in DESeq2 as it needs interger value
count_data<- round(count_data)
target_data<- prep_targets(target_matrix, delim = "c")

## Set significance cutoffs for p-values and log fold change
init_cutoffs(p_signif = 0.05, lfc_cut = 1.5)

## Initialize output mode
set_output_mode("both")

#### Preparing composite field for exp.design
## Creating a composite field for desired factors need to be analyzed, here, for example, "type" and "condition"
head(target_data)
target_data<- make_composite_field(c("condition", "type"))
head(target_data)

### Filtering was not utilized in this script

## Preparing DESeq2 object based on exp.design
dds<- prep_dds_from_data(count_input = count_data,
                         target_input = target_data,
                         experiment_design = ~ condition_type,
                         stabilization = "vst")

## Examining the data as a whole using Hierarchical Clustering
plot_euclid_dist(row_labels = "condition_type", filename = "Euclidian_distance.pdf", theme = 2, returnData = FALSE)

plot_dendro(filename = "Dendogram_plot.pdf", id_field = "targetID", groupBy = "condition" )
plot_dendro(filename = "Dendogram_plot.pdf", id_field = "targetID", groupBy = "type" )

## Multi-dimensional scaling analysis (MDS) (applicable for two factors)
plot_mds(filename = "MDS.pdf", color_var = "condition", shape_var = "type", theme = 1)

## Creating MDS Hull plot to identify outlier samples (this step is not needed if no outlier is observed in previous steps)
plot_mds_hulls(filename = "MDS_Hull.pdf",
               color_var = "condition",
               shape_var = "none",
               deOnly = FALSE,
               showLabel = TRUE,
               hullType = "solid",
               theme = 1)
# without sample labelling
plot_mds_hulls(filename = "MDS_Hull.pdf",
               color_var = "condition",
               shape_var = "type",
               deOnly = FALSE,
               showLabel = FALSE,
               hullType = "solid",
               theme = 1)

# viewing mds hull in another way
plot_mds_hulls(filename = "MDS_Hull_bytype.pdf",
               color_var = "type",
               shape_var = "condition",
               deOnly = FALSE,
               showLabel = FALSE,
               hullType = "solid",
               theme = 1)

## Visualization of Normalization 
plot_group_stats(filename = "Normalized_counts.pdf",
                 id_field = "targetID",
                 groupBy = "condition_type",
                 normalized = TRUE,
                 theme = 2)

#### Analyzing Differential Expression
dds<- DESeq(dds)

## Creating contrasts based on two factors using in this example: condition (Late and early as control) and type 
res_Late_case_typeA_vs_early_control<-results(dds,contrast = c("condition_type", "Late_A", "Early_A"))
res_Late_case_typeB_vs_early_control<-results(dds,contrast = c("condition_type", "Late_B", "Early_B"))

## Making a list of all the contrasts created 
result_list<- list(res_Late_case_typeA_vs_early_control,res_Late_case_typeB_vs_early_control)

## Combining differentially expressed genes across all contrast results
master_dataframe<- create_master_res(result_list, filename = "master_DE_list.txt", method = "union", lfc_filter =TRUE  )

## counting DE genes and visualization
de_counts(result_list, filename = "DE_counts.pdf",
          theme = 1)

## density plot against p-adjusted value
de_density_plot(result_list, filename = "aggregate_density.pdf",
                type = "pval",
                method = "union")

#### Examination of Differentially Expressed Genes
## Differences between Conditions
de_diverge_plot(result_list, filename = "DE_diverge.pdf", theme = 1, returnData = TRUE) ## not worked

de_boxplot(result_list, filename = "DE_boxplot.pdf", theme = 4) ## not worked


de_volcano(result_list, filename = "DE_volcano.pdf", theme = 1, strict_scale = TRUE)

## Differences Between Genes
de_heat(result_list, anno_columns = c("condition", "type"), filename = "upReg_heatmap.pdf", sort_choice = "max", numGenes = 25, theme = 2)

de_profile_plot(result_list, filename = "DE_profile_upReg_10.pdf",
                sort_choice = "max",
                numGenes = 25,
                theme = 1)

## To look at specific genes (did not work)
de_heat(result_list, anno_columns = c("condition", "type"),
        filename = "upReg_heatmap.pdf", sort_choice = "max",
        specific_genes = c("Zm00001d020704", "Zm00001d002649", "Zm00001d011461","Zm00001d031759", "Zm00001d043242"),
        cluster_contrasts = FALSE, theme = 2)


## To look at specific gene expression in boxplot
plot_gene(filename = "Zm00001d043242.pdf",
          gene_name = "Zm00001d043242", groupBy ="condition_type", theme = 4)

## Looking group of genes co-expression pattern
de_series(result_list, filename = "series_pattern.pdf", designVar = "condition_type",
          groupBy = "condition_type", method = "mean", numGroups = 4, writeData = TRUE, returnData = FALSE, theme = 1)
# with glm method
de_series(result_list, filename = "series_pattern.pdf", designVar = "condition_type",
          groupBy = "condition_type", method = "glm", writeData = TRUE, returnData = FALSE, theme = 1)
