## Loading necessary library
library(readxl)
library(dplyr)
library(writexl)
library(purrr)
## loading data
df<- read_excel("ad_for_structure.xlsx")

## replacing AA format to A/A format
df<- df %>%
  mutate(across(-`rs#`, ~ gsub("AA", "A/A", .))) %>% ##"~" is called anonymous function similar to function(); " ``" backticks are used to handle special characters in col name
  mutate(across(-`rs#`, ~ gsub("CC", "C/C", .))) %>%
                mutate(across(-`rs#`, ~gsub("GG", "G/G", .))) %>%
                              mutate(across(-`rs#`, ~gsub("TT", "T/T", .)))

## replacing A/A format to number format for structure software
df_left_replace<- df %>%
  mutate(across(-`rs#`, ~ gsub("/A", "", .))) %>% 
  mutate(across(-`rs#`, ~ gsub("/C", "", .))) %>%
  mutate(across(-`rs#`, ~gsub("/G", "", .))) %>%
  mutate(across(-`rs#`, ~gsub("/T", "", .)))

## replacing A/A format to number format for structure software
df_right_replace<- df %>%
  mutate(across(-`rs#`, ~ gsub("A/", "", .))) %>% 
  mutate(across(-`rs#`, ~ gsub("C/", "", .))) %>%
  mutate(across(-`rs#`, ~gsub("G/", "", .))) %>%
  mutate(across(-`rs#`, ~gsub("T/", "", .)))

## To find out if there are any unwanted letters present in the data, like "N"
unique_letter<- df_left_replace %>%
  select(-`rs#`) %>%
  map(unique) %>%
  unlist() %>%
  unique()
##
unique_letter_right<- df_right_replace %>%
  select(-`rs#`) %>%
  map(unique) %>%
  unlist() %>%
  unique()

## replacing letter format to number format for structure software
df_left_number<- df_left_replace %>%
  mutate(across(-`rs#`, ~ gsub("A", "1", .))) %>% 
  mutate(across(-`rs#`, ~ gsub("C", "2", .))) %>%
  mutate(across(-`rs#`, ~gsub("G", "3", .))) %>%
  mutate(across(-`rs#`, ~gsub("T", "4", .)))

## 
df_right_number<- df_right_replace %>%
  mutate(across(-`rs#`, ~ gsub("A", "1", .))) %>% 
  mutate(across(-`rs#`, ~ gsub("C", "2", .))) %>%
  mutate(across(-`rs#`, ~gsub("G", "3", .))) %>%
  mutate(across(-`rs#`, ~gsub("T", "4", .)))

## converting number to numeric as number is saved as text.
df_left_number<- df_left_number %>%
  mutate(across(-`rs#`, as.numeric))

df_right_number<- df_right_number %>%
  mutate(across(-`rs#`, as.numeric))

## saving the result in excel
write_xlsx(df_left_number, "number_format_left_letter_replace.xlsx")
        
write_xlsx(df_right_number, "number_format_right_letter_replace.xlsx")

## transposing final data for being ready for STRUCTURE
data<- read_excel("ad_data_ready_for_transpose.xlsx")

t_data<- t(data) ### transpose create matrix format

## converting matrix format to dataframe as my data is mixed type
t_data<- as.data.frame(t_data)

## saving the result in txt file
write.table(t_data, file = "ad_data_STRUCTURE.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
