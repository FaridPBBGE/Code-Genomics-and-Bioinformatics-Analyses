## Loading libraries
library(tidyverse)
library(readxl)
  graphics.off() 

# 1. Load and Transform
df <- read_excel("anox_aca1_for_analysis.xlsx")

df_long <- df %>%
  pivot_longer(
    cols = 3:74, 
    names_to = "SNP_ID", 
    values_to = "Allele"
  )

# 2. Setup Grid
snp_per_page <- 12
all_snps <- unique(df_long$SNP_ID)
snp_batch <- split(all_snps, ceiling(seq_along(all_snps) / snp_per_page))

# 3. Opening the PDF
pdf("ACA1_SNP_Score_result.pdf", width = 11, height = 8.5)

for (batch in snp_batch) {
  temp_df <- df_long %>% filter(SNP_ID %in% batch)
  
  # Ensure Score is numeric for the Y-axis (to show boxes)
  temp_df$Score <- as.numeric(temp_df$Score)
  
  # plotting
  plot <- ggplot(temp_df, aes(x = Allele, y = Score)) +
    theme_classic() +
    
    # BOXPLOT (Statistical summary)
    geom_boxplot(outlier.shape = NA, alpha = 0.5, fill="grey90", color = "black") +
    
    # JITTER (Colored by Phenotype Group)
    geom_jitter(aes(color = factor(Score)), width = 0.2, alpha = 0.5, size = 1.2) +
    
    # Grid
    facet_wrap(~SNP_ID, ncol = 3, scales = "free_x") +
    
    # Visual styling
    ##scale_fill_manual(values = c("A"="#E41A1C", "T" = "#377EB8", ##"C" = "#4DAF4A", "G" = "#984EA3")) +
    scale_color_manual(values = c("1"="#000000", "2"="cyan1", "3"="#D55E00")) +
    scale_y_continuous(breaks = c(1, 2, 3), limits = c(0.5, 3.5)) +
    
    labs(x = "Allele", y = "Anoxia Score (1-3)", color = "Phenotype Group") +
    
    # --- LEGEND ENHANCEMENT START ---
    guides(color = guide_legend(override.aes = list(size = 5))) + # Makes dots in legend bigger
    theme(
      strip.text = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold"),
      
      # Bolding Allele letters
      axis.text.x = element_text(face = "bold", color = "black"),
      
      # Legend Text and Title Sizes
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"), # Legend Title size
      legend.text = element_text(size = 12),                # Legend Labels size
      legend.key.size = unit(1.2, "cm")                     # Spacing between legend items
    )
  # --- LEGEND ENHANCEMENT END ---
  
  print(plot)
}

dev.off()


#=======================================================================#
## Creating box plots in a single page

#=======================================================================#
## Loading libraries
library(tidyverse)
library(readxl)
graphics.off() 

# 1. Load and Transform
df <- read_excel("cax1_for_analysis.xlsx")

df_long <- df %>%
  pivot_longer(
    cols = 3:74, 
    names_to = "SNP_ID", 
    values_to = "Allele"
  ) %>%
  mutate(Score = as.numeric(Score))

## Open TIFF device with high resolution (300 DPI)
tiff(
  filename = "Cax1_SNP_score_1page_result.tiff",
  width = 18,        # Width in inches (adjust if needed)
  height = 14,       # Height in inches
  units = "in", 
  res = 300,         # Publication quality DPI
  compression = "lzw" # Reduces output file size significantly
)

## Create single-page plot with tight spacing
plot <- ggplot(df_long, aes(x = Allele, y = Score)) +
  theme_classic() +
  
  # BOXPLOT (Statistical summary)
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill="grey90", color = "black", lwd = 0.3) +
  
  # JITTER (Colored by Phenotype Group)
  geom_jitter(aes(color = factor(Score)), width = 0.2, alpha = 0.6, size = 0.8) +
  
  # Grid (9 columns x 8 rows = 72 subplots)
  facet_wrap(~SNP_ID, ncol = 9, scales = "free_x") +
  
  # Visual styling
  scale_color_manual(values = c("1" = "#000000", "2" = "cyan1", "3" = "#D55E00")) +
  scale_y_continuous(breaks = c(1, 2, 3), limits = c(0.5, 3.5)) +
  
  labs(x = "Allele", y = "Anoxia Score (1-3)", color = "Phenotype Group") +

# Legend and margin styling
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(
    # --- MARGINS & SPACING ---
    plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "in"), # 0.1 inch edge margins
    panel.spacing = unit(0.15, "lines"),                                   # Tight gap between subplots
    
    # --- TEXT SCALING FOR 72 FACETS ---
    strip.text = element_text(face = "bold", size = 7),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", color = "black", size = 6.5),
    axis.text.y = element_text(size = 6),
    
    # --- LEGEND ---
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.5, "cm")
  )
print(plot)

dev.off()


#====Getting Statistics result for these SNPs from GWAS output file=======

#=========================================================================

## Loading libraries
library(dplyr)
library(readxl)
library(writexl)

## load data
master_file<- read.csv("GAPIT.Association.GWAS_Results.MLMM.Anox.csv")
snp_list<- read_excel("anox_aca1_snp_score_data.xlsx", sheet = "SNP")

## Subsetting inf. from the master file based on targeted SNP list
subset_SNP_file<- left_join(snp_list, master_file, by = "SNP")

## Save the subset file in excel
write_xlsx(subset_SNP_file, "ACA1_SNP_statistics_inf.xlsx")
