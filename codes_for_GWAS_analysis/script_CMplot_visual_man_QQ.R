## Installing CMplot package
install.packages("CMplot")

## Loading the package
library(CMplot)

## Loading necessary packages
library(readxl)
library(readr)
library(dplyr)

## Loading data
data<- read_xlsx("AD_Anox_manhat_noGEMMA.xlsx")

## Creating SNP-density plot
CMplot(data, plot.type = "d", bin.size =1e4, chr.den.col = c("darkgreen", "yellow", "red"), file = "jpg", file.name ="Arabidopsis", dpi = 300,  file.output = TRUE, verbose = TRUE, width = 9, height = 6)


## Creating stacked manhattan plot for various models together Genome-wide
CMplot(data, plot.type = "m", multracks =TRUE, col = c("grey60", "#4197d8"), threshold =1e-4, threshold.lty =2, threshold.lwd = 2, threshold.col = "red", amplify = TRUE, signal.col = "red",signal.cex = 1, file = "jpg", file.name = "AD_Ca_multi_manhattn_Plot", file.output = TRUE, verbose = TRUE, width =16)


## Creating stacked manhattan plot for various models together chromosome-wise
CMplot(data, plot.type = "m", multracks =TRUE, col = c("#4197d8"), amplify = FALSE, file = "tiff", file.name = "Anox_Chr5_manhattn_Plot_highlightsnp", file.output = TRUE, verbose = TRUE, width =16)

## Creating top 100 SNPs based on top p-value of model per chromosome and highlight them in the plot
top_cmlm<- data %>% 
  group_by(Chromosome) %>%
  arrange(CMLM) %>%
  slice(1:100) %>%
  ungroup() %>%
  pull(SNP)

top_bl<- data %>% 
  group_by(Chromosome) %>%
  arrange(BLINK) %>%
  slice(1:100) %>%
  ungroup() %>%
  pull(SNP)

#top_gm<- data %>% 
#  arrange(GEMMA) %>%
 # slice(1:100) %>%
  #pull(SNP)

top_ml<- data %>% 
  group_by(Chromosome) %>%
  arrange(MLMM) %>%
  slice(1:100) %>%
  ungroup() %>%
  pull(SNP)

# Store the resutls in a list
SNPs<- list(
  top_cmlm,
top_bl,
#top_gm,
top_ml
  )


## Creating stacked manhattan plot for various models together chromosome-wise with highlighted SNPs (good picture one--use setting of this for other figures)
CMplot(data, plot.type = "m", multracks =TRUE, col = c("#4197d8"), amplify = FALSE, width =16, highlight=SNPs, highlight.col="red", axis.cex = 1.5, legend.cex =2.0, file = "jpg", file.name = "Anox_Chr5_man_Plot_highlightsnp", file.output = TRUE, verbose = TRUE,)


## Creating stacked manhattan plot for various models together with highlighting top SNPs Genome-wide
CMplot(data, plot.type = "m", multracks = TRUE, col = c("grey60", "#4197d8"), amplify =FALSE, width =18, highlight = SNPs, highlight.col = "red", axis.cex = 1.8, legend.cex =2, file = "jpg", file.name = "Anox_50ksnp_100hlightsnp", file.output = TRUE, verbose = TRUE)

!!!!!~~~~~~~Creating Zoomed Manhattan plot-Specific genomic region~~~~~!!!!!

## Loading Whole genome data
gen_data<- read_excel("AD_Anox_manhat_noGEMMA.xlsx")

## Defining the speficic genomic region to be explored
zoom_start<- 6539423
zoom_end<- 6549012

## subsetting zoomed genomic region data from the main data file
zoomed_data<- gen_data[gen_data$Chromosome==2 &
                       gen_data$Position >= zoom_start &
                       gen_data$Position <= zoom_end, ]

## Setting up highlighting SNPs and labeling (here, SNPs located in selective genes will be colored in green and signficant SNPs that are absent in gene will be colored as red)

# data reading for getting name of targeted SNPs
target_snps_file<- read_excel("list_snp_highlight_zoom_man_plot.xlsx", sheet = "chr2-g2")

# Manually indexing green and red SNP position
green_idx<-c(1)
red_idx<- c(2, 3, 4, 5)

# Extract the SNP names
green_snps<- target_snps_file$SNP[green_idx]
red_snps<- target_snps_file$SNP[red_idx]

# Extract the gene IDs for labels (green SNPs only)
green_snp_labels<- target_snps_file$GeneID[green_idx]

## Combining these two colored SNPs list for highlighting in the plot
all_highlighted_snps<- c(green_snps, red_snps)

## Setting up labeling (I will label those green SNPs with gene names, where those SNPs are located)

# Setting up for automatic labeling and coloring for next steps
num_green<- length(green_snps)
num_red<- length(red_snps)

label_snps<- c(green_snp_labels, rep("", num_red)) # follow the same order of color snps list joined in all_highlighted_snps variable

## Creating the manhattan plot for targeted genomic region
CMplot(zoomed_data[,c(1,2,3,5)],
       plot.type = "m",
       col = c("grey60"),
       cex = 1.4,
       
	   #Highlight settings
       highlight=all_highlighted_snps,
	     highlight.cex = 1.4,
       highlight.col=c(rep("green", num_green), rep("red", num_red)),
	     highlight.pch = 19,
       
	   #Text label settings
       highlight.text = label_snps,
       highlight.text.col = "black",
       highlight.text.cex = 1.2,
       
    # Layout settings
       amplify = FALSE,
       main = "Chr2: 6.54-6.55 Mb",
       ylab = expression(-log[10](italic(p)) * "-BLINK"),
       file.output = TRUE,
       file = "jpg",
       file.name = "chr2_g2_blink_zoomed_5kb",
       verbose = TRUE,
       width =12,
       height=6)



#============== Creating manhattan plot (chromosome-wise with highlight and labeling SNPs) ===================
  
## Loading required packages
library(CMplot)
library(readxl)
library(dplyr)

## Loading data 
data <- read_xlsx("anox_Chr_3.xlsx")

#========================================================================
# STEP 1: Extract Top 100 Highlights Per Model
# ========================================================================
# --- CMLM ---
# Sort by lowest p-value and extract top 100 for highlights
top_cmlm<- data%>% arrange(CMLM) %>% slice(1:100) %>% pull(SNP)

# --- BLINK ---
top_bl<- data%>% arrange(BLINK) %>%  slice(1:100) %>% pull(SNP)

# --- MLMM ---
top_ml <- data%>% arrange(MLMM) %>%  slice(1:100) %>% pull(SNP)


## This keeps only the global top 100 SNPs that are physically inside in Chr5 file
top_cmlm_chr5 <- top_cmlm_genome[top_cmlm_genome %in% data$SNP]
top_bl_chr5   <- top_bl_genome[top_bl_genome %in% data$SNP]
top_ml_chr5   <- top_ml_genome[top_ml_genome %in% data$SNP]


#======= Extract Top 100 Highlights Per model per Chr and Save the result in excel for downstream analysis ========#
# Load library
library(openxlsx)

data_Chr<- read_xlsx("anox_Chr_1.xlsx")

#Extract Top 100 Highlights Per model
top_cmlm_100<- data_Chr %>% 
  arrange(CMLM) %>% 
  slice(1:100) %>% 
  select(SNP, Chromosome, Position) %>% 
  mutate(Model = "CMLM")

top_bl_100<- data_Chr %>% 
  arrange(BLINK) %>% 
  slice(1:100) %>% 
  select(SNP, Chromosome, Position) %>% 
  mutate(Model = "BLINK")

top_mlmm_100<- data_Chr %>% 
  arrange(MLMM) %>% 
  slice(1:100) %>% 
  select(SNP, Chromosome, Position) %>% 
  mutate(Model = "MLMM")

# Combine the result and Format for downstream analaysis (common model will be together with "comma")
# Stacking all three tables on top of each other
combined_snps<- bind_rows(top_cmlm_100, top_bl_100, top_mlmm_100)

# Group by SNP name and combing/collapse overlapped model names
result_excel<- combined_snps %>%
  group_by(SNP,Chromosome, Position) %>%
  summarise(MOdel_name= paste(Model, collapse = ";"), .groups = "drop") %>%
  arrange(Chromosome, Position)

# Save the result in excel
write.xlsx(result_excel, file = "anox_Chr1_top100_SNP.xlsx")

#============== Plotting =========================

#=================================================

#-----------  CMLM ----------------------
Chr<-3

CMplot(data %>% select(SNP, Chromosome, Position, CMLM), 
           plot.type = "m", 
           multracks = FALSE,                 
           col = c("#4197d8"), 
           amplify = FALSE, 
           width = 16, 
           highlight = top_cmlm,       
           highlight.col = "red", 
           highlight.text = top_cmlm,  
           highlight.text.col = "black",      
           highlight.text.cex = 0.4,          
           axis.cex = 1.5, 
           legend.cex = 2.0, 
           file = "pdf", 
           file.name = paste0("Anox_Chr", Chr, "_CMLM_0.4"), 
           file.output = TRUE, 
           verbose = TRUE)
  
#----------- BLINK ---------------------

CMplot(data %>% select(SNP, Chromosome, Position, BLINK), 
       plot.type = "m", 
       multracks = FALSE,                 
       col = c("#4197d8"), 
       amplify = FALSE, 
       width = 16, 
       highlight = top_bl,       
       highlight.col = "red", 
       highlight.text = top_bl,  
       highlight.text.col = "black",      
       highlight.text.cex = 0.4,          
       axis.cex = 1.5, 
       legend.cex = 2.0, 
       file = "pdf", 
       file.name = paste0("Anox_Chr", Chr, "_BLINK_0.4"), 
       file.output = TRUE, 
       verbose = TRUE)


#----------- MLMM ----------------------
 
CMplot(data %>% select(SNP, Chromosome, Position, MLMM), 
       plot.type = "m", 
       multracks = FALSE,                 
       col = c("#4197d8"), 
       amplify = FALSE, 
       width = 16, 
       highlight = top_ml,       
       highlight.col = "red", 
       highlight.text = top_ml,  
       highlight.text.col = "black",      
       highlight.text.cex = 0.4,          
       axis.cex = 1.5, 
       legend.cex = 2.0, 
       file = "pdf", 
       file.name = paste0("Anox_Chr", Chr, "_mlmm_0.4"), 
       file.output = TRUE, 
       verbose = TRUE)


#=============================== Done ====================================#
  
  

##===================QQ Plot visualization===========================##

## Multiple QQ plots all together in one plot(final)
data<- read_csv("AD__allsnp_Anox_nogemma_manhat_QQ_view.csv")

CMplot(data, plot.type = "q", col= c("#4197d8", "darkgoldenrod1","darkkhaki"), multraits =TRUE, ylab.pos = 2,box = FALSE, legend.cex =2, axis.cex = 1.5, file = "tiff", file.name = "AD_Anox_QQ_allsnp", dpi = 300, file.output =TRUE, verbose =TRUE, ylim = c(0,10), width = 8, height = 6, conf.int=FALSE)


## Multiple QQ plots separately but in the same figure
data<- read_xlsx("AD_Ca_manhat_QQ_view.xlsx")

CMplot(data, plot.type = "q", col= c("dodgerblue1", "darkgoldenrod1","aquamarine","darkkhaki"), multracks =TRUE, ylab.pos = 2, threshold.col ="red", threshold.lty = 2, conf.int = TRUE, box = FALSE, axis.cex = 1, file = "jpg", file.name = "AD_Ca_QQ1",threshold = 1e-4, dpi = 300, file.output =TRUE, verbose =TRUE, ylim = c(0,8), width = 6, height = 6)


## Single QQ plot using for loop for assigning a specific color for downstream analysis

models<- c("CMLM", "BLINK", "MLMM")
colors<- c("#4197d8", "darkgoldenrod1","darkkhaki")

## looping through each model

for (m in seq_along(models)) {
 
   CMplot(data %>% select(SNP, Chromosome, Position, all_of(models[m])),
         plot.type = "q", 
         col= colors[m],
         ylab.pos = 2,
         box = FALSE,
         legend.cex =1.8,
         axis.cex = 1.5,
         file = "tiff",
         file.name = paste0("AD_Anox_QQ_", models[m]),
         dpi = 300,
         file.output =TRUE,
         verbose =TRUE,
         ylim = c(0,10),
         width = 8,
         height = 6,
         conf.int=FALSE)
}


