## Loading library
library(dplyr)
library(readr)
library(readxl)
library(writexl)
library(tidyr)
library(ggplot2)
library(purrr)

## reading/uploading the files into R
F1<- read_excel("TPM_result__visual.xlsx")
F2<- read_excel("TPM_TPMCalculator_V4_maize_Hi-A.xlsx")

## renaming 1st column genes names in F2 to match gene's name in F1
F2$Gene<- gsub("\\.RefGen_V4$","",F2$Gene)

## Sub-setting F2 file based on Gene list in F1 file using dplyr
F2_subset<- F2 %>% 
  filter(Gene %in% F1$Gene)

## Renaming the columns names in F2 file to match with the cols names in F1 file
colnames(F2_subset)[-1]<- colnames(F1)[-1]

## Combining the two files
combined_F<- bind_rows(F1,F2_subset)

## Adding a new variable based on F1 and F2 information
combined_F$Tools<- c(rep("Salmon", nrow(F1)), rep("TPMCalc", nrow(F2_subset)))

## Converting "NA" strings to actual NA values
combined_F[combined_F== "NA"]<- NA

## Converting all columns to numeric except Gene and Tools cols
combined_F<- combined_F %>%
  mutate(across(-c(Gene, Tools), ~ as.numeric(.)))


## Visualization with ggplot
# To visualize in ggplot, data should be in long format; not in wide format.
long_F_df<- combined_F %>%
  pivot_longer(cols = -c("Gene","Tools"), names_to = "Sample", values_to = "TPM")

## Checking missing values
missing_value<- long_F_df %>% 
  filter(is.na(TPM))


## Counting the value based on tool type for each sample
count_df<- long_F_df %>%
  group_by(Sample, Tools) %>%
  summarise(Count= sum(!is.na(TPM)))

## Pivoting data to wider format for correlation
wide_F_df <- long_F_df %>%
  pivot_wider(names_from = Tools, values_from = TPM)

## Filtering out rows with missing values
wide_F_df <- wide_F_df %>%
  filter(!is.na(Salmon) & !is.na(TPMCalc))

## Calculating correlatino
correlations <- wide_F_df %>%
  group_by(Sample) %>%
  summarize(correlation = cor(Salmon, TPMCalc, use = "complete.obs"))

## merging correlation result to wide format data
wide_F_df<- wide_F_df %>%
  left_join(correlations, by = "Sample")

## missing value counting based on Tools type for each sample
na_count_df<- long_F_df %>%
  group_by(Sample, Tools) %>%
  summarise(NA_Count= sum(is.na(TPM)))


## Visualization of the counting result
ggplot(count_df, aes(x= Tools, y=Count, fill = Tools)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = Count), vjust= -0.5, size =3) +
  facet_wrap(~ Sample, scales = "free_y") +
  theme_minimal() +
  labs(x="Tools", y="Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
## Visualization for correlation
ggplot(wide_F_df, aes(x = Salmon, y = TPMCalc)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Sample, scales = "free") +
  labs(title = "Correlation of TPM values between Salmon and TPMCalc",
       x = "TPM (Salmon)",
       y = "TPM (TPMCalc)") +
  theme_minimal()
