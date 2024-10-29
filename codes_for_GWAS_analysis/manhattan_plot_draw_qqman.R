## opening the qqman
library(qqman)
library(readxl)

### Q-Q plot
Q_Plot<- read_excel("C:/Users/Farid/Documents/R/R_output/RS_qq_manhattan_plot.xlsx")

qq(Q_Plot$P, pch=19, col="blue", cex=2)

### Manhattan plot
df_manhat<- read_excel("manhattan_plot_Ca_gemma.xlsx")
manhattan(df_manhat,col=c("palegreen2","skyblue"),ylim=c(0,8),cex=2, genomewideline = -log10(7.297e-05),suggestiveline = FALSE,highlight = NULL, logp = TRUE, annotatePval = NULL)

manhattan(df_manhat,col=c("palegreen2","skyblue"),ylim=c(0,8),cex=2, genomewideline = -log10(8.108e-06),suggestiveline = FALSE,highlight = NULL, logp = TRUE, annotatePval = NULL)

qq(df_manhat$P, pch=19, col="blue",cex=2)
