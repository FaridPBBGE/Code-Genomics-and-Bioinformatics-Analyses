####loading software
library(mrMLM)
### Loading data
library(readxl)

### codes for multi-locus GWAS

mrMLM(fileGen = "C:/Users/Farid/Documents/R/R_output/Multilocus_models_SNP.txt", filePhe = "C:/Users/Farid/Documents/R/R_output/K_multilocus.txt", fileKin = NULL, filePS =  "C:/Users/Farid/Documents/R/R_output/multi_Q_matrix.txt", PopStrType = "Q", fileCov = NULL, Genformat = "Hmp", method = c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),Likelihood =  "REML", trait = 1, SearchRadius = 20, CriLOD = 3, SelectVariable = 50, Bootstrap = FALSE, DrawPlot = TRUE, Plotformat = "jpeg", dir = "C:/Users/Farid/Documents/R/R_output")

#### csv file format
mrMLM(fileGen = "C:/Users/Farid/Documents/R/R_output/Minerals_models_SNP.csv", filePhe = "C:/Users/Farid/Documents/R/R_output/K_multilocus.csv", fileKin = NULL, filePS =  "C:/Users/Farid/Documents/R/R_output/multi_Q_matrix.csv", PopStrType = "Q", fileCov = NULL, Genformat = "Hmp", method = c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),Likelihood =  "REML", trait = 1, SearchRadius = 20, CriLOD = 3, SelectVariable = 50, Bootstrap = FALSE, DrawPlot = TRUE, Plotformat = "jpeg", dir = "C:/Users/Farid/Documents/R/R_output")

### To draw plot using results file from mrMLM 
MultiManhattan(ResultIntermediate ="C:/Users/Farid/Documents/R/R_output/1_intermediate result.csv",ResultFinal = "C:/Users/Farid/Documents/R/R_output/1_Final result.csv", mar = c(2.5, 2.5, 0.7, 2.5),PlotFormat = "tiff","jpeg")

library("mrMLM.GUI")
mrMLM.GUI()
