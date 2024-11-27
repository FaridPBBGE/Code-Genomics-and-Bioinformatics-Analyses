## Installing CMplot package
install.packages("CMplot")

## Loading the package
library(CMplot)

## Loading necessary packages
library(readxl)
library(readr)
library(dplyr)

## Loading data
  data<- read_xlsx("AD_Zn_manhat_QQ_view.xlsx")

## Creating SNP-density plot
CMplot(data, plot.type = "d", bin.size =1e4, chr.den.col = c("darkgreen", "yellow", "red"), file = "jpg", file.name ="Arabidopsis", dpi = 300,  file.output = TRUE, verbose = TRUE, width = 9, height = 6)

## Creating stacked manhattan plot for various models together
CMplot(data, plot.type = "m", multracks =TRUE, col = c("grey60", "#4197d8"), threshold =1e-4, threshold.lty =2, threshold.lwd = 2, threshold.col = "red", amplify = TRUE, signal.col = "red",signal.cex = 1, file = "jpg", file.name = "AD_Ca_multi_manhattn_Plot", file.output = TRUE, verbose = TRUE, width =16)

## Creating top 100 SNPs based on top p-value of model and highlight them in the plot
top_cmlm<- data %>% 
  arrange(CMLM) %>%
  slice(1:1000) %>%
  pull(SNP)

top_bl<- data %>% 
  arrange(BLINK) %>%
  slice(1:1000) %>%
  pull(SNP)

top_gm<- data %>% 
  arrange(GEMMA) %>%
  slice(1:1000) %>%
  pull(SNP)

top_ml<- data %>% 
  arrange(MLMM) %>%
  slice(1:1000) %>%
  pull(SNP)

# Store the resutls in a list
SNPs<- list(
  top_cmlm,
top_bl,
top_gm,
top_ml
  )

## Creating stacked manhattan plot for various models together with highlighting top SNPs
CMplot(data, plot.type = "m", multracks = TRUE, col = c("grey60", "#4197d8"), amplify = TRUE, signal.col = "red",signal.cex = 1, file = "jpg", file.name = "AD_Zn", file.output = TRUE, verbose = TRUE, width =18, highlight = SNPs)

## Multiple QQ plots all together in one plot
data<- read_xlsx("AD_Anox_manhat_QQ_view.xlsx")

CMplot(data, plot.type = "q", col= c("dodgerblue1", "darkgoldenrod1","aquamarine","darkkhaki"), multraits =TRUE, threshold = 1e-4, ylab.pos = 2, threshold.col ="red", threshold.lty = 2, conf.int = TRUE, box = FALSE, axis.cex = 1, file = "jpg", file.name = "AD_Anox_QQ", dpi = 300, file.output =TRUE, verbose =TRUE, ylim = c(0,8), width = 7, height = 6)

## Multiple QQ plots separately but in the same figure
data<- read_xlsx("AD_Ca_manhat_QQ_view.xlsx")

CMplot(data, plot.type = "q", col= c("dodgerblue1", "darkgoldenrod1","aquamarine","darkkhaki"), multracks =TRUE, ylab.pos = 2, threshold.col ="red", threshold.lty = 2, conf.int = TRUE, box = FALSE, axis.cex = 1, file = "jpg", file.name = "AD_Ca_QQ1",threshold = 1e-4, dpi = 300, file.output =TRUE, verbose =TRUE, ylim = c(0,8), width = 6, height = 6)
