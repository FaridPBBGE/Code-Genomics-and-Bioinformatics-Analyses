## Necessary libraries
library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
#################### Visualizing Fst result ######################
# Read data
fst_data <- read_delim("100K_Anox_intol__vs_Anox_tol.windowed.weir.fst", delim = "\t")

# Create midpoint position
fst_data <- fst_data %>%
  mutate(MID_POS = (BIN_START + BIN_END) / 2)

### Manhattan Plot — Chromosome-wise (Faceted)
# 95th percentile threshold
fst_threshold <- quantile(fst_data$WEIGHTED_FST, 0.90, na.rm = TRUE)

# Manhattan plot faceted by chromosome
ggplot(fst_data, aes(x = MID_POS, y = WEIGHTED_FST, color = as.factor(CHROM))) +
  geom_point(size = 1) +
  geom_hline(yintercept = fst_threshold, linetype = "dashed", color = "red") +
  facet_wrap(~ CHROM, scales = "free_x") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Position within Chromosome (bp)", y = "Weighted Fst",
       title = "Fst (Chromosome-wise)") +
  theme(
    panel.background = element_blank(),        # no gray background
    panel.grid.major = element_blank(),        # no major grid lines
    panel.grid.minor = element_blank(),        # no minor grid lines
    axis.line = element_line(color = "black"), # add axis lines
    strip.background = element_blank(),        # remove facet label background
    legend.position = "none"
  )

### Cumulative Genome-wide Position with ggplot2
# Compute chromosome offsets
chr_offsets <- fst_data %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(BIN_END)) %>%
  mutate(chr_start = cumsum(lag(chr_len, default = 0)))

# Merge with data and compute cumulative position
fst_data <- fst_data %>%
  left_join(chr_offsets, by = "CHROM") %>%
  mutate(CUM_POS = MID_POS + chr_start)

# Manhattan plot with cumulative positions
ggplot(fst_data, aes(x = CUM_POS, y = WEIGHTED_FST, color = as.factor(CHROM))) +
  geom_point(size = 1) +
  geom_hline(yintercept = fst_threshold, linetype = "dashed", color = "red") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Cumulative Genomic Position (bp)", y = "Weighted Fst",
       title = "Fst Manhattan Plot (Whole Genome - Cumulative)") +
  theme(
    panel.background = element_blank(),         # white panel background
    panel.grid.major = element_blank(),         # remove major grid lines
    panel.grid.minor = element_blank(),         # remove minor grid lines
    axis.line = element_line(color = "black"),  # clean axis lines
    legend.position = "none"
  )


### running in a loop
# List all .fst files
fst_files <- list.files(pattern = "\\.fst$")

# Function to process and plot each Fst file
plot_fst_manhattan <- function(file) {
  # Read the data
  fst_data <- read_delim(file, delim = "\t", show_col_types = FALSE)
  
  # Calculate midpoint and cumulative position
  fst_data <- fst_data %>%
    mutate(MID_POS = (BIN_START + BIN_END) / 2)
  
  chr_offsets <- fst_data %>%
    group_by(CHROM) %>%
    summarise(chr_len = max(BIN_END)) %>%
    mutate(chr_start = cumsum(lag(chr_len, default = 0)))
  
  fst_data <- fst_data %>%
    left_join(chr_offsets, by = "CHROM") %>%
    mutate(CUM_POS = MID_POS + chr_start)
  
  # Calculate 95th percentile threshold
  fst_threshold <- quantile(fst_data$WEIGHTED_FST, 0.90, na.rm = TRUE)
  
  # Create plot
  p <- ggplot(fst_data, aes(x = CUM_POS, y = WEIGHTED_FST, color = as.factor(CHROM))) +
    geom_point(size = 1) +
    geom_hline(yintercept = fst_threshold, linetype = "dashed", color = "red") +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Cumulative Genomic Position (bp)", y = "Weighted Fst",
         title = paste("Genome-wide Plot -", file)) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none"
    )
  
  # Save plot as JPEG
  ggsave(filename = paste0(tools::file_path_sans_ext(file), "_manhattan.jpeg"),
         plot = p,
         width = 10, height = 4, dpi = 300)
}

# Apply the function to all files
lapply(fst_files, plot_fst_manhattan)


### Smoothed Line Plot Function
fst_files <- list.files(pattern = "\\.fst$")

# Smoothed Line Plot Function
plot_fst_smoothed <- function(file) {
  # Read data
  fst_data <- read_delim(file, delim = "\t", show_col_types = FALSE)
  
  # Calculate midpoint
  fst_data <- fst_data %>%
    mutate(MID_POS = (BIN_START + BIN_END) / 2)
  
  # Compute cumulative genome-wide position
  chr_offsets <- fst_data %>%
    group_by(CHROM) %>%
    summarise(chr_len = max(BIN_END)) %>%
    mutate(chr_start = cumsum(lag(chr_len, default = 0)))
  
  fst_data <- fst_data %>%
    left_join(chr_offsets, by = "CHROM") %>%
    mutate(CUM_POS = MID_POS + chr_start)
  
  # 95th percentile threshold for Fst
  fst_threshold <- quantile(fst_data$WEIGHTED_FST, 0.90, na.rm = TRUE)
  
  # Smoothed plot
  p <- ggplot(fst_data, aes(x = CUM_POS, y = WEIGHTED_FST, color = as.factor(CHROM))) +
    geom_smooth(method = "loess", se = FALSE, span = 0.3, size = 0.7) +  # Smoothing line
    geom_hline(yintercept = fst_threshold, linetype = "dashed", color = "red") +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Cumulative Genomic Position (bp)", y = "Weighted Fst",
         title = paste("Smoothed Fst Plot -", file)) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none"
    )
  
  # Save JPEG
  ggsave(filename = paste0(tools::file_path_sans_ext(file), ".jpeg"),
         plot = p,
         width = 10, height = 4, dpi = 300)
}

# Apply to all files
lapply(fst_files, plot_fst_smoothed)


################# for pi result #############
library(tools) 
pi_files <- list.files(pattern = "\\.pi$")

plot_pi_cumulative <- function(file) {
  # Read data
  pi_data <- read_delim(file, delim = "\t", show_col_types = FALSE)
  
  # Midpoint of each window
  pi_data <- pi_data %>%
    mutate(MID_POS = (BIN_START + BIN_END) / 2)
  
  # Compute chromosome offsets for cumulative position
  chr_offsets <- pi_data %>%
    group_by(CHROM) %>%
    summarise(chr_len = max(BIN_END)) %>%
    mutate(chr_start = cumsum(lag(chr_len, default = 0)))
  
  pi_data <- pi_data %>%
    left_join(chr_offsets, by = "CHROM") %>%
    mutate(CUM_POS = MID_POS + chr_start)
  
  # Top 10% threshold
  pi_thresh <- quantile(pi_data$PI, 0.90, na.rm = TRUE)
  
  # Create plot
  p <- ggplot(pi_data, aes(x = CUM_POS, y = PI, color = as.factor(CHROM))) +
    geom_smooth(method = "loess", span = 0.3, se = FALSE, size = 0.7) +
    geom_hline(yintercept = pi_thresh, linetype = "dashed", color = "red") +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Cumulative Genome Position (bp)",
         y = "Nucleotide Diversity (π)",
         title = paste("Genome-wide π -", file)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  # Save JPEG
  filename_out <- paste0(file_path_sans_ext(file), "_pi_plot.jpeg")
  ggsave(filename = filename_out, plot = p, width = 10, height = 4, dpi = 300)
}

# Apply the function to all pi files
lapply(pi_files, plot_pi_cumulative)


################### Visualizing xpclr result  ###############################
# Loading library
library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tools)

# Getting all xpclr files
xpclr_files<- list.files(pattern = "^anox_Chr_")

# Creating function to plot smoothed line grpah for xpclr with two dashed lines (one for mean and another for top 5% value)

plot_smoothed_xpclr<- function(file) {
# reading data
data<- read_delim(file, delim = "\t", show_col_types = FALSE)

# calculating midpoint
data<- data %>% mutate(mid_pos=(start + stop)/2)

# calculating 95th percentile or top 5% value and mean of the data
  xpclr_thresh<- quantile(data$xpclr, 0.95,na.rm= TRUE)
  xpclr_mean<- mean(data$xpclr, na.rm=TRUE)
        
# plotting
grpah<- ggplot(data, aes(x= mid_pos, y= xpclr)) +
geom_line(color="steelblue", alpha=0.5) +
geom_hline(yintercept = xpclr_thresh, linetype= "dotted", color="red") +
geom_hline(yintercept = xpclr_mean, linetype= "dashed", color="blue") +
labs(
 x= "Genomic Position",
y= "XPCLR"
) +
theme_classic() +
theme(
panel.grid = element_blank(),
 axis.line = element_line(color = "black")
 )
# Saving the plot
result<- paste0(file, "_xpclr_smoothed.jpeg")
ggsave(result, plot =grpah, width = 10, height = 4, dpi= 300)
}
## Running the function for all the files
lapply(xpclr_files, plot_smoothed_xpclr)


################### Combined two plots in the same figure using package ###########

##### To combine two images in the same figure, file's names should be similar except file extension. I am creating images using loop, so, I am using "plot_smoothed_xpclr" and
"plot_smoothed_fst" functions with saving images as object, instead of saving as images #####

##### For Fst result using loop

### Creating function to plot smoothed line grpah for Fst files with two dashed lines (one for mean and another for top 5% value)
plot_smoothed_fst<- function (file) {
# reading data
data_fst<- read_delim(file, delim= "\t", show_col_types=FALSE)

# calculating midpoint position
data_fst <- data_fst %>%
  mutate(MID_POS = (BIN_START + BIN_END)/2)
		
# calculating 95th percentile or top 5% value and mean of the data
fst_thresh<- quantile(data_fst$WEIGHTED_FST, 0.95, na.rm = TRUE)
fst_mean<- mean(data_fst$WEIGHTED_FST, na.rm=TRUE)

# plotting
plot_fst<- ggplot(data_fst, aes(x= MID_POS, y= WEIGHTED_FST)) +
geom_line(color="steelblue", alpha=0.5) +
geom_hline(yintercept = fst_thresh, linetype= "dotted", color="red") +
geom_hline(yintercept = fst_mean, linetype= "dashed", color="blue") +
labs(
 x= "Genomic Position",
y= "Weighted_Fst"
) +
theme_classic() +
theme(
panel.grid = element_blank(),
 axis.line = element_line(color = "black")
 )
return(plot_fst)

}

## Plotting all Fst files as an object

#Getting all Fst files
Fst_files<- list.files(pattern = "\\.fst$")

# plotting all files together
fst_plots<- setNames(
  lapply(Fst_files, plot_smoothed_fst),
  tools::file_path_sans_ext(Fst_files)
)

######## For xpclr result using loop

# Creating function to plot smoothed line grpah for xpclr with two dashed lines (one for mean and another for top 5% value)

plot_smoothed_xpclr<- function(file) {
# reading data
data<- read_delim(file, delim = "\t", show_col_types = FALSE)

# calculating midpoint
data<- data %>% mutate(mid_pos=(start + stop)/2)

# calculating 95th percentile or top 5% value and mean of the data
  xpclr_thresh<- quantile(data$xpclr, 0.95,na.rm= TRUE)
  xpclr_mean<- mean(data$xpclr, na.rm=TRUE)
        
# plotting
plot_xpclr<- ggplot(data, aes(x= mid_pos, y= xpclr)) +
geom_line(color="steelblue", alpha=0.5) +
geom_hline(yintercept = xpclr_thresh, linetype= "dotted", color="red") +
geom_hline(yintercept = xpclr_mean, linetype= "dashed", color="blue") +
labs(
 x= "Genomic Position",
y= "XPCLR"
) +
theme_classic() +
theme(
panel.grid = element_blank(),
 axis.line = element_line(color = "black")
 )
return(plot_xpclr)
}

## Plotting all xpclr files as an object

# Getting all xpclr files
xpclr_files<- list.files(pattern = "^Zn_Chr_\\d+$")

# plotting all files together
xpclr_plots<-setNames(
  lapply(xpclr_files, plot_smoothed_xpclr),
  xpclr_files
  )

#### Creating combined plots using cowplot package
library(cowplot)

common_names<- intersect(names(fst_plots), names(xpclr_plots))

for (name in common_names) {
  combined_plot<- plot_grid(
    fst_plots[[name]],
    xpclr_plots[[name]],
    ncol = 1,
    align = "v",
    rel_heights = c(1,1)
  )
  ggsave(
    filename = paste0(name, "_fst_xpclr.jpeg"),
    plot = combined_plot,
    width = 12, height = 6, dpi = 300
  )
}


############ Extracting top 5% Fst and xp-clr value ##################

#### Using function to extract top 5% value for multiple files
library(readr)
library(dplyr)
library(writexl)
library(tools)

## Creating function
subset_5percent<- function(data) {

# getting base file name
base_name_file<-file_path_sans_ext(basename(data)) 

# reading data
df<- read_delim(file = data, delim = "\t", show_col_types = FALSE)
  
# filtering data based on criteria
  subset_df<- df %>%
    filter(xpclr >= quantile(xpclr,probs=0.95, na.rm=TRUE))

# saving the result
  file_name<- paste0(base_name_file, "_top5%_xpclr.xlsx")
  write_xlsx(subset_df,file_name)
}

##  applying the function for all files

# getting files, where function will be applied
list_files<- list.files(path = ".", pattern = "_Chr_\\d+$")

# applying the function
lapply(list_files,subset_5percent)

