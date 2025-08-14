## installing tximport
BiocManager::install("tximport")
## loading packages
library(tximport)
library(dplyr)
library(ggplot2)
library(stringr)
library(writexl)
## getting files
c_files<- list.files(path =".", pattern = ".sf", full.names = TRUE, recursive = TRUE)
## renaming files with accession number
names(c_files)<- stringr::str_split(c_files,pattern = "/", simplify = TRUE)[,2] %>% ##[] will select which column I need to select for replace. The column will be output of str_split code
stringr::str_replace("_quant","") ## removing quant 


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
  dplyr::filter(type== "mRNA")%>% ### for gtf file, mRNA will be transcript
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
write_xlsx(tpm_gene_df,"TPM_result_Salmon_maize_V4.xlsx")
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
#######################0########################################################
################################################################################

## DEseq2 for multi-factor experiment 

library(DESeq2)

meta_data_mf<-data.frame(condition=c("Early", "Early", "Early", "Early", "Late", "Late", "Late", "Late"), type=c("A","A","B","B","A","A","B","B"))

## naming object file's row name according to the columns of counts
rownames(meta_data_mf)<-colnames(t2gene_test$counts)

## Converting characters to factor
meta_data_mf$type<- factor(meta_data_mf$type)

dea_mf<- DESeqDataSetFromTximport(t2gene_test,colData= meta_data_mf, design = ~ condition)
dea_mf$condition<- relevel(dea$condition, ref = "Early")
dea_1<-DESeq(dea_mf)
result_1<- results(dea_1, contrast = c("condition","Late", "Early" ))
result_1

## adding another factor (later stage) to the model
deamf_add<- dea_1
levels(deamf_add$type)

## running DeSeq in another way#1st
design(deamf_add)<- formula(~ type + condition)
deamf_add<- DESeq(deamf_add)

## running DeSeq with multi-factor#2nd
dea_mf_alt<- DESeqDataSetFromTximport(t2gene_test,colData= meta_data_mf, design = ~ type + condition)
dea_mf_alt$condition<- relevel(dea$condition, ref = "Early")
dea_1_alt<-DESeq(dea_mf_alt)
result_1_alt<- results(dea_1_alt, contrast = c("condition","Late", "Early" ))
result_1_alt

## getting results
result_mf<- results(deamf_add)
head(result_mf)

## Looking 2nd variable independently likewise 1st variable
result_mf_type<- results(deamf_add, contrast = c("type", "A", "B"))
head(result_mf_type)


## To know overall comparison test done by Deseq2
resultsNames(deamf_add)


## Count data transformation
vst_d<-vst(deamf_add, blind = FALSE)
rlog_d<- rlog(deamf_add, blind = FALSE)
head(assay(vst_d),5)

## Transformation's effect on the variance and visualization
normal_ctf<- normTransform(deamf_add) ## using log2(n+1)

## For visualization these count transformed data
BiocManager::install("vsn")
library("vsn")
meanSdPlot(assay(normal_ctf))
meanSdPlot(assay(vst_d))
meanSdPlot(assay(rlog_d))

### Checking data quality by sample clustering and visualization
## Heatmap for the count matrix
library(pheatmap)
select_top_exgene<- order(rowMeans(counts(deamf_add, normalized=TRUE)),decreasing = TRUE)[1:20]
df_dea<- as.data.frame(colData(deamf_add)[,c("condition","type")])

## For normal log transformation
pheatmap(assay(normal_ctf)[select_top_exgene,],cluster_rows = FALSE, show_rownames = FALSE,cluster_cols = FALSE, annotation_col = df_dea)

## For vsd transformation
pheatmap(assay(vst_d)[select_top_exgene,],cluster_rows = FALSE, show_rownames = FALSE,cluster_cols = FALSE, annotation_col = df_dea)

## For rlog transformation
pheatmap(assay(rlog_d)[select_top_exgene,],cluster_rows = FALSE, show_rownames = FALSE,cluster_cols = FALSE, annotation_col = df_dea)

## Heatmap for sample to sample distances
Dist_samples<- dist(t(assay(vst_d)))
library(RColorBrewer)

## Converting distance object to matrix
Dist_samples_matrix<- as.matrix(Dist_samples)
rownames(Dist_samples_matrix)<- paste(vst_d$condition, vst_d$type, sep = "-")
colnames(Dist_samples_matrix)<- NULL
color<- colorRampPalette(rev(brewer.pal(7,"Greens")))(255)
pheatmap(Dist_samples_matrix,
         clustering_distance_rows = Dist_samples,
         clustering_distance_cols = Dist_samples,
         col=color)

## PCA analysis of the samples
plotPCA(vst_d,intgroup=c("condition","type"))

## PCA Visualization in ggplot
d_plotPCA<- plotPCA(vst_d,intgroup=c("condition","type"), returnData=TRUE, ntop=2000)
variance_percent<- round(100 * attr(d_plotPCA, "percentVar"))
ggplot(d_plotPCA, aes(PC1,PC2, color=condition, shape = type)) +
  geom_point(size=3)+
  xlab(paste0("PC1: ",variance_percent[1],"% variance"))+
  ylab(paste0("PC2: ",variance_percent[2],"% variance")) +
  coord_fixed()

##### Some variation in the standard workflow (for in case situation) 

## To test Wald test
dds_wald_test<- estimateSizeFactors(deamf_add)
dds_wald_test<- estimateDispersions(deamf_add)
dds_wald_test<- nbinomWaldTest(deamf_add)
results(dds_wald_test)

## More topics on shrinkage estimators
### Three methods are available: apeglm, ashr,and normal method
## for ageglm method
result_aptT<-lfcShrink(deamf_add, coef = 3, type = "apeglm", lfcThreshold = 1)
plotMA(result_aptT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="red", lwd=2)

## for ashr method
result_ashT<-lfcShrink(deamf_add, coef = 2, type = "ashr", lfcThreshold = 1)
plotMA(result_ashT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)

## Inspecting outlier counts
par(mar=c(8,5,2,2))
boxplot(log10(assays(deamf_add)[["cooks"]]), range=0, las=2)

## Dispersion plot and gene-wise dispersion estimate
plotDispEsts(deamf_add)

## visualization of independent filtering result done by DESeq2
metadata(result_mf)$alpha
metadata(result_mf)$filterThreshold
plot(metadata(result_mf)$filterNumRej,
     type="b", ylab="Rejections Number",
     xlab="Quantiles of filter")
lines(metadata(result_mf)$lo.fit, col="red")
abline(v=metadata(result_mf)$filterTheta)

## log2 fold change testing with different hypothesis

res_greaterAbs<-results(deamf_add,lfcThreshold = 0.5, altHypothesis = "greaterAbs") ## for two-tailed test

res_lessAbs<-results(deamf_add,lfcThreshold = 0.5, altHypothesis = "lessAbs") ## for maximum of the upper and lower tests

res_greater<-results(deamf_add,lfcThreshold = 0.5, altHypothesis = "greater") ## for  testing > than lfcthreshold value

res_less<-results(deamf_add,lfcThreshold = 0.5, altHypothesis = "less") ## for testing < lfcthreshold value

## visualization of alternative tests results
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim<-c(-3,3)
drawlines<- function() abline(h=c(-0.5,0.5), col="dodgerblue", lwd=2)
plotMA(res_greaterAbs, ylim=ylim); drawlines()
plotMA(res_lessAbs, ylim=ylim); drawlines()
plotMA(res_greater, ylim=ylim); drawlines()
plotMA(res_less, ylim=ylim); drawlines()

## Looking into calculated values done by DESeq2
mcols(deamf_add, use.names = TRUE)[1:5,1:5]
## To know the description of the columns 
mcols(mcols(deamf_add), use.names = TRUE)[1:5,]
## To see all the names of the columns 
names(mcols(deamf_add)) 
## To see mean value for each gene and samples
head(assays(deamf_add)[["mu"]])
## To see cook's distance value for each gene and samples
head(assays(deamf_add)[["cooks"]])
## To check dispersion
head(dispersions(deamf_add)) or
head(mcols(deamf_add)$dispersion)
## To check size factors
sizeFactors(deamf_add)
## coefficient check
head(coef(deamf_add))
## To know general information about prior used for lfc shrinkage
priorInfo(file_object)

## Count's outlier inspection
Wald_stat<-result_mf$stat
max_cook<- apply(assays(deamf_add)[["cooks"]],1, max)
index_wald_stat<- !is.na(Wald_stat)
plot(rank(Wald_stat[index_wald_stat]),max_cook[index_wald_stat], xlab="Rank of Wald statistic",ylab="Max. Cook's distance per gene",
     ylim=c(0,5), cex=0.4, col=rgb(0,0,1,.3))
g<-2
l<-ncol(deamf_add)
abline(h=qf(.95,g, l-g), col="red")

## p-value observation
plot(result_mf$baseMean+1, -log10(result_mf$pvalue),log="x", xlab="Mean of normalized counts",
     ylab= expression(-log[10](pvalue)),ylim=c(0,30),cex=0.4, col=rgb(0,0,1,0.3))

## p-value distribution
filtered_pvalue<- result_mf$baseMean > metadata(result_mf)$filterThreshold
histgram1<- hist(result_mf$pvalue[!filtered_pvalue],breaks = 0:50/50, plot = FALSE)
histgram2<- hist(result_mf$pvalue[filtered_pvalue],breaks = 0:50/50, plot = FALSE)
color_pattern<- c('do not pass'="Khaki", 'pass'="powderblue")

barplot(height = rbind(histgram1$counts,histgram2$counts),beside = FALSE,col = color_pattern, space = 0, main = "p-value distribution",ylab = "frequency")
text(x=c(0, length(histgram1$counts)),y=0, label=paste(c(0,1)),adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(color_pattern), legend = rev(names(color_pattern)))

## Counts plot
plotCounts(deamf_add, gene = which.min(result_mf$padj),intgroup = c("condition", "type"))

## Plotting counts for two factors of variables (applicable time-series exp as well, time will be in x-axis)
count_plot2fac<- plotCounts(deamf_add, gene = which.min(result_mf$padj), intgroup = c("condition","type"), returnData = TRUE)

## Plotting 
ggplot(count_plot2fac,
       aes(x=condition, y=count, color=type, group = type)) + geom_point() +
         stat_summary(fun = mean, geom = "line") +
         scale_y_log10()

#### Clustering significant genes based on variables or factors of study
parameters<- coef(deamf_add)
colnames(parameters)
# Plotting significant genes
sig_genes<- head(order(result_mf$padj),20)
matrix_genes<- parameters[sig_genes,-c(1)]
thr_limit<- 3
matrix_genes[matrix_genes< -thr_limit]<- -thr_limit
matrix_genes[matrix_genes>thr_limit]<- thr_limit

library(pheatmap)
pheatmap(matrix_genes, breaks = seq(from= -thr_limit, to=thr_limit, length=101),
         cluster_cols = FALSE)

#### Plotting differential expressed genes onto the genome ####

## For genome features, gff or gtf file will be used as Granges required for karyoplot package
library(rtracklayer)
library(dplyr)

file_path<- file.path(".","Zmays_RefGen_V4.gene.gff3")## giving file path
f<-rtracklayer::import(file_path)
head(f)

## converting gff3 file to dataframe to extract the target inf.
f_df<- as.data.frame(f)
head(f_df)
## getting gene name from this df file, so that it can be used with DEseq2 result obj., where same gene name format used
gene_list <- f_df %>%
  filter(type == "gene") %>%
  select(-source, -score, -phase, -ID, -pacid, -longest, -Parent)

## To remove gene column
gene_list <- gene_list %>%
  select(-type) %>%
  rename(gene_id=Name)

head(gene_list)

## Converting the data to GRanges object
library(GenomicRanges)

gr_gene_db<- GRanges(seqnames = gene_list$seqnames, ranges = IRanges(start = gene_list$start, end = gene_list$end), strand = gene_list$strand)

mcols(gr_gene_db)<- data.frame(gene_id= gene_list$gene_id)

head(gr_gene_db)

## Adding some colums from DESeq2 result object
names(gr_gene_db)<- gr_gene_db$gene_id

## Checking if common gene IDs between DeSeq2 result obj and newly created GRanges obj.
common_genes<- intersect(names(gr_gene_db), rownames(result_mf))

## Adding cols from DeSeq2 obj to GRanges obj for only common genes
mcols(gr_gene_db)<- result_mf[names(gr_gene_db),c("log2FoldChange", "stat", "pvalue", "padj")]
head(gr_gene_db)

## To create custom genome for plotting in karyotype package,  df of gff file can be used to get necessary inf. for creating custom genome for any species using the respective species's gtf or gff file
my_genes<- f_df %>%
  filter(type=="gene") %>%
  select(seqnames, start, end) %>%
  rename(chr=seqnames)
head(my_genes)

## Creating custom genome based on chr length
library(readxl)
chr_length<- read_excel("Zmays_V4_chr_length.xlsx", col_names = TRUE)
custom_genome<- toGRanges(data.frame(
  chr=chr_length$Chromosome,
  start=1,
  end=chr_length$Size
))

## Creating custom cytobands
custom_cytoband<- toGRanges(data.frame(
  chr=my_genes$chr,
  start=my_genes$start,
  end=my_genes$end,
  name=paste0(my_genes$chr,"_gene", seq_len(nrow(my_genes))),
  gieStain=ifelse((my_genes$end - my_genes$start) >20000, "gpos100",                 ifelse((my_genes$end - my_genes$start) > 10000, "gpos75",
           ifelse((my_genes$end - my_genes$start) > 5000, "gpos50", "gpos25")))
  ))

## Alternative way to create cytoband for less dense
custom_cytoband_alt<- do.call(rbind, lapply(seq_len(nrow(my_genes)), function(i) {
  chr_name<- my_genes$chr[i]
  chr_start<- my_genes$start[i]
  chr_end<- my_genes$end[i]
  num_bands<- 10
  
  # Dividing the genomic region into bands
  data.frame(
    chr= chr_name,
    start= seq(chr_start, chr_end, length.out= num_bands +1)[-num_bands - 1],
    end=seq(chr_start, chr_end, length.out= num_bands + 1)[-1],
    name= paste0(chr_name,"_band", seq_len(num_bands)),
    gieStain=sample(c("gneg", "gpos25", "gpos50", "gpos75","gpos100"),num_bands, replace = TRUE, prob = c(0.4, 0.2, 0.2, 0.1, 0.1))
    )
  
}))

## Converting to GRanges
custom_cytoband_alt<- toGRanges(custom_cytoband_alt)



## Plotting using Karyotype and using custom genome
library(karyoploteR)
kp<- plotKaryotype(genome = custom_genome, cytobands = custom_cytoband)

## Plotting top 10 genes names that differentially expressed
gene_ordered<- gr_gene_db[order(gr_gene_db$padj, na.last = TRUE),]
kp<- kpPlotMarkers(kp, gene_ordered[1:10], labels = names(gene_ordered[1:10]), text.orientation = "horizontal")

## Filtering out NA in padj column and convert them to -log10 to visualize them easily
filtered_gr_gene_db<- gr_gene_db[!is.na(gr_gene_db$padj)]
log.padj<- -log10(filtered_gr_gene_db$padj)
mcols(filtered_gr_gene_db)$log.padj<- log.padj
head(filtered_gr_gene_db)

## Plotting the only significant genes below padj < 0.05 value
sig_genes_karyo<- filtered_gr_gene_db[filtered_gr_gene_db$padj < 0.05]

## visualizing log2FoldChange over the chromosome
#To see the range value for log2FoldChange
range(sig_genes_karyo$log2FoldChange)

# As lfc value has above and below 0 (positive and negative value), it is better to make make Y-axis symmetric around 0 for visual clarity.
lfc.ymax<- ceiling(max(abs(range(sig_genes_karyo$log2FoldChange))))
lfc.ymin<- -lfc.ymax

# To visualize differential expressed genes with increasing/decreasing dot size 
plot<- getDefaultPlotParams(plot.type = 2)
plot$data2height<- 75
plot$ideogramheight<- 14

kp<- plotKaryotype(genome = custom_genome, cytobands = custom_cytoband, plot.type = 2, plot.params = plot)

# To add Title
kpAddMainTitle(kp, main = "Differential Gene Expression")


## Data for plot gene expression
dot_size<- sqrt(sig_genes_karyo$log.padj)/6

# Adding color to differentiate between over and underexpressed genes
col_over<- "#FFBD07AA"
col_under<- "#00A6EDAA"
sig_col<- ifelse(sig_genes_karyo$log2FoldChange < 0, col_under,col_over)

kpPoints(kp, data = sig_genes_karyo, cex = dot_size, y=sig_genes_karyo$log2FoldChange, ymax = lfc.ymax, ymin = lfc.ymin, r1=points_top, col=sig_col)

# If no need to add axis and labels, no need to run kpaxis and kpAddLabels (optional)
kpAxis(kp, ymax = lfc.ymax, ymin = lfc.ymin, cex=0.8, r1=points_top)
kpAddLabels(kp, labels = "logfc", srt=90, pos = 1, label.margin = 0.04, ymax= lfc.ymax, ymin=lfc.ymin, r1=points_top)

# Adding top genes names in the plot
top_genes<- gene_ordered[1:10]

# To adjust the labelling
points_top<- 0.8

kpPlotMarkers(kp, top_genes, labels = names(top_genes), text.orientation = "horizontal",label.color = "blue",line.color = "#777777", cex=0.6, r0=points_top)

gene_midpoint<- start(top_genes) + (end(top_genes) - start(top_genes))/2 
(## gene.midpoint = center point of each gene (start + half the gene's length)

kpSegments(kp, chr= as.character(seqnames(top_genes)), x0=gene_midpoint, x1=gene_midpoint, y0=top_genes$log2FoldChange, y1=lfc.ymax, ymax = lfc.ymax, ymin = lfc.ymin, r1=points_top, col="#777777")

## Data for gene density
# Adding gene density below the ideograms
kp<- kpPlotDensity(kp, data = gr_gene_db, window.size = 50e4, data.panel = 2)
