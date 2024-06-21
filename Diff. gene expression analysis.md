# Differential gene expression analysis using DESeq2 tool
**updating**
~~~
## installing tximport
BiocManager::install("tximport")
## loading packages
library(tximport)
library(dplyr)
library(ggplot2)
library(stringr)
## getting files
c_files<- list.files(path =".", pattern = ".sf", full.names = TRUE, recursive = TRUE)
## renaming files with accession number
names(c_files)<- stringr::str_split(c_files,pattern = "/", simplify = TRUE)[,2] %>% ##[] will select which column I need to select for replace. The column will be output of str_split code
stringr::str_replace("_quant","")


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
## getting transcript and gene ID from this df file, so that it can be used for salmon count files
t2gene<- f_gff3_df%>%
  dplyr::filter(type== "mRNA")%>%
  select(ID,Parent)%>%
  mutate(ID=str_replace(ID,pattern= "\\.RefGen_V4$",""),
         Parent=str_replace(Parent,pattern="\\.RefGen_V4$",""))

## importing t2gene file into tximport
t2gene_test<- tximport(c_files,type = "salmon",tx2gene = t2gene)

## installing DESeq2
BiocManager::install("DESeq2")
## loading DESeq2
library(DESeq2)

meta_data<-data.frame(condition=c("Early", "Early", "Early", "Early", "Late", "Late", "Late", "Late"))

rownames(meta_data)<-colnames(t2gene_test$counts)

## DESeq2 workflow
dea<- DESeqDataSetFromTximport(t2gene_test,colData= meta_data, design = ~ condition)
dea$condition<- relevel(dea$condition, ref = "Early")
dea<-DESeq(dea)
result<- results(dea, contrast = c("condition","Late", "Early" ))
result
## Log fold change shrinkage for visualization and ranking
resultsNames(dea)
BiocManager::install("apeglm")
library(apeglm)

result_lfc_sh<- lfcShrink(dea, coef = "condition_Late_vs_Early", type = "apeglm")

result_lfc_sh
## Ordering p-values and adjusted p-values
result_ord<- result[order(result$pvalue),]

## Summarize the result
summary(result)

## To know how many adjusted p-values are < 0.1
sum(result$padj < 0.1, na.rm=TRUE)
### To know details function on results, use ?results
## To change alpha value=0.05, default is 0.1
result_05<-results(dea,alpha = 0.05)

## Summarize the result with new alpha value
summary(result_05)

#### To know how many adjusted p-values are < 0.05
sum(result$padj < 0.05, na.rm=TRUE)

## To visualize the result-MA-Plot
# first trying with normal log2 fold change result
plotMA(result, ylim=c(-2,2))
# now trying with shrinkage log fold change
plotMA(result_lfc_sh, ylim=c(-2,2))

## To detect row number of genes by clicking and save the result
det_r_gene<- identify(result$baseMean, result$log2FoldChange)
rownames(result)[det_r_gene]

## To calculate alternative shrinkage method other than apeglm
result_lfc_nor<- lfcShrink(dea, coef =2, type = "normal")
result_lfc_ashr<- lfcShrink(dea, coef =2, type = "ashr")

## Visualize the all methods together
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim<-c(1,1e5); ylim<- c(-3, 3)
plotMA(result_lfc_sh, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(result_lfc_nor, xlim=xlim, ylim=ylim, main="normal")
plotMA(result_lfc_ashr, xlim=xlim, ylim=ylim, main="ashr")

## Counting reads per group plotting (plotting with default)
plotCounts(dea, gene = which.min(result$padj), intgroup = "condition")

## To plot counting read customize way using ggplot
c_plot<- plotCounts(dea, gene = which.min(result$padj), intgroup = "condition", returnData = TRUE)

ggplot(c_plot,aes(x=condition,y=count))+
  geom_point(position = position_jitter(w=0.1, h=0))+
  scale_y_log10(breaks=c(100,300,600, 900))

## To know more informaiton about result object got before
mcols(result)$description
~~~

