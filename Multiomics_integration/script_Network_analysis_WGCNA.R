library(tidyverse)
library(writexl)
library(readxl)
ori_data<- read_excel("Gene_expression-matrix_ZnEx-Root.xlsx")
### view(column_names)
column_names
# Calculate means for specific groups
data_means <- ori_data %>%
  mutate(
   # Leaf_cont_0h_mean = rowMeans(select(., starts_with("Leaf_cont_0h_"))),
   # Leaf_cont_2d_mean = rowMeans(select(., starts_with("Leaf_cont_2d_"))),
    #Leaf_FeLim_2d_mean = rowMeans(select(., starts_with("Leaf_FeLim_2d_"))),
    #Leaf_cont_4d_mean = rowMeans(select(., starts_with("Leaf_cont_4d_"))),
    # Leaf_FeLim_4d_mean = rowMeans(select(., starts_with("Leaf_FeLim_4d_"))),
    # Leaf_cont_7d_mean = rowMeans(select(., starts_with("Leaf_cont_7d_"))),
    # Leaf_FeLim_7d_mean = rowMeans(select(., starts_with("Leaf_FeLim_7d_"))),
    # Leaf_cont_14d_mean = rowMeans(select(., starts_with("Leaf_cont_14d_"))),
    #Leaf_FeLim_14d_mean = rowMeans(select(., starts_with("Leaf_FeLim_14d_"))),
    # Leaf_cont_21d_mean = rowMeans(select(., starts_with("Leaf_cont_21d_"))),
    # Leaf_FeLim_21d_mean = rowMeans(select(., starts_with("Leaf_FeLim_21d_")))
   # Root_cont_0h_mean = rowMeans(select(., starts_with("Root_cont_0h_"))),
    #Root_cont_1h_mean = rowMeans(select(., starts_with("Root_cont_1h_"))),
    #Root_FeLim_1h_mean = rowMeans(select(., starts_with("Root_FeLim_1h_"))),
    #Root_cont_2d_mean = rowMeans(select(., starts_with("Root_cont_2d_"))),
    #Root_FeLim_2d_mean = rowMeans(select(., starts_with("Root_FeLim_2d_"))),
   # Root_cont_4d_mean = rowMeans(select(., starts_with("Root_cont_4d_"))),
    #Root_FeLim_4d_mean = rowMeans(select(., starts_with("Root_FeLim_4d_"))),
    #Root_cont_7d_mean = rowMeans(select(., starts_with("Root_cont_7d_"))),
    #Root_FeLim_7d_mean = rowMeans(select(., starts_with("Root_FeLim_7d_"))),
    #Root_cont_14d_mean = rowMeans(select(., starts_with("Root_cont_14d_"))),
    #Root_FeLim_14d_mean = rowMeans(select(., starts_with("Root_FeLim_14d_"))),
   # Root_cont_21d_mean = rowMeans(select(., starts_with("Root_cont_21d_"))),
   # Root_FeLim_21d_mean = rowMeans(select(., starts_with("Root_FeLim_21d_")))
	
    Leaf_cont_0h_mean = rowMeans(select(., starts_with("Leaf_cont_0h_"))),
   # Leaf_cont_2d_mean = rowMeans(select(., starts_with("Leaf_cont_2d_"))),
    #Leaf_ZnEx_2d_mean = rowMeans(select(., starts_with("Leaf_ZnEx_2d_"))),
    #Leaf_cont_4d_mean = rowMeans(select(., starts_with("Leaf_cont_4d_"))),
    # Leaf_ZnEx_4d_mean = rowMeans(select(., starts_with("Leaf_ZnEx_4d_"))),
    # Leaf_cont_7d_mean = rowMeans(select(., starts_with("Leaf_cont_7d_"))),
    # Leaf_ZnEx_7d_mean = rowMeans(select(., starts_with("Leaf_ZnEx_7d_"))),
    # Leaf_cont_14d_mean = rowMeans(select(., starts_with("Leaf_cont_14d_"))),
    #Leaf_ZnEx_14d_mean = rowMeans(select(., starts_with("Leaf_ZnEx_14d_"))),
    # Leaf_cont_21d_mean = rowMeans(select(., starts_with("Leaf_cont_21d_"))),
    # Leaf_ZnEx_21d_mean = rowMeans(select(., starts_with("Leaf_ZnEx_21d_")))
    Root_cont_0h_mean = rowMeans(select(., starts_with("Root_cont_0h_"))),
    Root_cont_1h_mean = rowMeans(select(., starts_with("Root_cont_1h_"))),
    Root_ZnEx_1h_mean = rowMeans(select(., starts_with("Root_ZnEx_1h_"))),
    Root_cont_2d_mean = rowMeans(select(., starts_with("Root_cont_2d_"))),
    Root_ZnEx_2d_mean = rowMeans(select(., starts_with("Root_ZnEx_2d_"))),
    Root_cont_4d_mean = rowMeans(select(., starts_with("Root_cont_4d_"))),
    Root_ZnEx_4d_mean = rowMeans(select(., starts_with("Root_ZnEx_4d_"))),
    Root_cont_7d_mean = rowMeans(select(., starts_with("Root_cont_7d_"))),
    Root_ZnEx_7d_mean = rowMeans(select(., starts_with("Root_ZnEx_7d_"))),
    Root_cont_14d_mean = rowMeans(select(., starts_with("Root_cont_14d_"))),
    Root_ZnEx_14d_mean = rowMeans(select(., starts_with("Root_ZnEx_14d_"))),
    Root_cont_21d_mean = rowMeans(select(., starts_with("Root_cont_21d_"))),
    Root_ZnEx_21d_mean = rowMeans(select(., starts_with("Root_ZnEx_21d_"))))
  

# Select the gene_id and the newly created mean columns
data_selected <- data_means %>%
  select(Gene_Id, ends_with("_mean"))

# Save the data frame to an Excel file
write_xlsx(data_selected, "mean_FeLim-Root_gene_matrix.xlsx")

## Loading data
df_data<- read_excel("mean_ZnEx-leaf_gene_matrix.xlsx")
## Renaming the header of the first column
names(df_data)[names(df_data)=="Gene_Id"]<-"Treatment"
write_xlsx(df_data, "mean_ZnEx_leaf_gene_matrix.xlsx")

### WGCNA loading
library(WGCNA)

# To transpose the data for converting it into matrix for loading data into WGCNA
df<- read_excel("mean_FeLim-leaf_gene_matrix.xlsx")


df_mat<- as.matrix(df[,-1]) ## To remove first column header;otherwise, it creates problem down analysis
row.names(df_mat)<- df$Treatment
t_df_mat<- t(df_mat)
t_df_mat[1:5,1:5]
## installing some dependencies packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"))

## selection soft threshold level
powers=c(c(1:10),seq(from=12, to=20, by=2))

soft_threshold<- pickSoftThreshold(t_df_mat,
                                   powerVector = powers,
                                   verbose = 5
                                   )

## visualize the soft_threshold result
par(mfrow=c(1,2));
cex1=0.9
plot(soft_threshold$fitIndices[,1],
     -sign(soft_threshold$fitIndices[,3])* soft_threshold$fitIndices[,2],
     xlab = "Soft Threshold (Power)",
     ylab = "Scale Free Topology Model Fit,R^2",
     main = paste("Scale Independence")
)
text(soft_threshold$fitIndices[,1],
     -sign(soft_threshold$fitIndices[,3])* soft_threshold$fitIndices[,2],
     labels = powers, cex = cex1, col = "blue", pos = 1
     )
abline(h=0.80, col="red")
plot(soft_threshold$fitIndices[,1],
     soft_threshold$fitIndices[,5],
     xlab = "Soft Threshold (Power)",
     ylab = "Mean Connectivity",
     main = paste("Mean Connectivity")
     )
text(soft_threshold$fitIndices[,1],
     soft_threshold$fitIndices[,5],
     labels = powers, cex = cex1, col = "blue",pos=1
)
abline(v=10, col="red")

## building network
picked_power=10
temp_cor<- cor
cor<- WGCNA::cor
network<- blockwiseModules(t_df_mat,
                           ## Adjacency function#
                           power = picked_power,
                           networkType = "signed",
                           maxBlockSize = 20000,
                           #Topological overlap options
                           TOMType = "signed",
                           #Basic Tree and Block Options#
                           deepSplit = 2,
                           pamRespectsDendro = F,
                           minModuleSize = 30,
                           #Module Adjustments#
                           mergeCutHeight = 0.25,
                           #TOM == Archive the run results in TOM file (saves time)
                           saveTOMs = T,
                           saveTOMFileBase ="FeLim_root",
                           #Output Options
                           numericLabels = T,
                           verbose = 3,
                           randomSeed=12345
                           )
cor<- temp_cor
library(readr)
readr::write_rds(network,
                 file = file.path("Fe_lim_treatment", "FeLim_root_treatment.RDS"))
## In case analysis cannot be done by a single day, network result can be save into RDS file and later can be opened by the following code.
net_treatment<-readRDS("Fe_lim_treatment/FeLim_leaf_treatment.RDS")
names(net_treatment)

## module review
mergedColors<- labels2colors(net_treatment$colors)
plotDendroAndColors(
  net_treatment$dendrograms[[1]],
  mergedColors[net_treatment$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  hang=0.05,
  addGuide = TRUE,
  guideHang = 0.05,
  cex.colorLabels =1.2,
  )

table_format<-net_treatment$colors[net_treatment$blockGenes[[1]]]
write.table(table_format, "genes_group_leaf")

my_table<- table(table_format)
my_table
my_table_df<-as.data.frame(my_table)
write_xlsx(my_table_df, "summary_genes_groups_leaf.xlsx")

## Listing module and corresponding genes list
list_modules<- data.frame(
  gene_ID= names(net_treatment$colors),
  colors=labels2colors(net_treatment$colors)
)

write_xlsx(list_modules, "List_genes_modules_root.xlsx")

## Getting module Eigengenes per cluster
eigen_module<- moduleEigengenes(t_df_mat, mergedColors)$eigengenes

## Sorting modules so similar modules are next to each other
eigen_module<- orderMEs(eigen_module)
module_sorted<- names(eigen_module)%>% gsub("ME","", .)

## Adding treatment names
eigen_module$treatment<- row.names(eigen_module)

## process data
v_eigen_module<- eigen_module%>%
  pivot_longer(-treatment) %>%
  mutate(
    name=gsub("ME","",name),
    name=factor(name,levels=module_sorted)
  )
## visulization_correlation between modules and traits
v_eigen_module%>% ggplot(.,aes(x=treatment,y=name,fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit=c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90, face="bold", size = 10), axis.text.y = element_text(face = "bold", size = 12)) +
  labs(title = "Module-treatment relationship", y="Modules", fill="corr")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  theme(axis.title.x = element_text(size = 16, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"))

## Making expression profiles of modules
row.names(list_modules)<-list_modules$gene_ID

df_module<- data.frame(df_mat)%>%
  mutate(
    gene_ID=row.names(.)
  )%>%
  pivot_longer(-gene_ID)%>%
  mutate(
    module= list_modules[gene_ID,]$colors
  )

## visualizing module exp.profiles
df_module%>% ggplot(., aes(x=name,y=value, group=gene_ID)) +
    geom_line(aes(color=module),
              alpha=0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90,face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)
  ) +
  facet_grid(rows = vars(module)) +
  labs(title="Modulewise gene expression profile",x="Stages of Treatment", y= "gene expression (TPM)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
  theme(axis.title.x = element_text(size = 16, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"))
 
### Creating TOMfile

TOM<- TOMsimilarityFromExpr(t_df_mat,power = picked_power)

row.names(TOM)=row.names(df_mat)
colnames(TOM)=row.names(df_mat)

edge_list<- data.frame(TOM)%>%
  mutate(
    gene1=row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2=name, correlation=value) %>%
  unique()%>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1=list_modules[gene1,]$colors,
    module2=list_modules[gene2,]$colors
  )

write_delim(edge_list,file = "fe_lim_edge_list.tsv", delim ="\t")
