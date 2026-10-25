## Loading required library
library(xml2)
library(writexl)

## Reading raw content as UTF-8 text, ignoring the declared encoding
raw_text <- readLines("panther_no_correction_analysis.xml", encoding = "UTF-8")

## Fixing the encoding declaration in the XML header
raw_text <- sub('encoding="UTF-16"', 'encoding="UTF-8"', raw_text, fixed = TRUE)

## Parse from the corrected text
doc <- read_xml(paste(raw_text, collapse = "\n"))

## Extracting all results blocks from the xml file
results<- xml_find_all(doc, ".//result")

## Initialize vectors
gene_ids<-c()
go_ids<- c()
go_labels<- c()
pvalues<- c()

## loop through each results
for (res in results) {
    go_id <- xml_text(xml_find_first(res, "./term/id"))
    go_label <- xml_text(xml_find_first(res, "./term/label"))
    pval<- xml_text(xml_find_first(res, ".//pValue"))
    #Each result block can have multiple mapped_id entries
    mapped_id<- xml_find_all(res, ".//mapped_id")
    #looping each mapped id for creating row for individual id
    for (map in mapped_id) {
      gene_ids<- c(gene_ids, xml_text(map))
      go_ids<- c(go_ids, go_id)
      go_labels<- c(go_labels, go_label)
      pvalues<- c(pvalues, pval)
    }
}

## Creating a dataframe to save the output
df<- data.frame(
  GeneID= gene_ids,
  GOterm= go_ids,
  Label= go_labels,
  pValue= pvalues,  stringsAsFactors = FALSE
)

## Saving as excel file
write_xlsx(df, "panther_gene_go_no_correction.xlsx")

~~~~~~~~~Merging this newly creating file with idmapping file, where Arabidopsis tair ID are converted to Uniprotkb ID~~~~~~~~~~~~~~~~~~

## loading libraries
library(dplyr)
library(readxl)
library(writexl)

## reading panther file
panther<- read.csv("pather_gene_go_no_correction.csv", stringsAsFactors = FALSE)

## reading Uniprot ID conversion file for merging wiht panther file
id<- read_excel("idmapping_2026_04_28.xlsx")

## keeping only necessary files from id file before merging
id<- id %>% select(
  "From", 
  "Entry", 
  "Reviewed",
  "Entry Name", 
  "Protein names",
  "Gene Names"
)

## rename the first column to match panther file
colnames(id)[1]<- "GeneID"

## merging two files
merged_file<- panther %>%
  left_join(id, by="GeneID")

## saving the merged_file
write_xlsx(merged_file, "merged_panther_uniprotid.xlsx")

## Extracting stringIDs that don't have input Tair gene ID, using mapping file to find out which stringid in interaction file does not have tair id using dplyr

## loading required libraries
library(dplyr)
library(readxl)
library(writexl)

## read interaction file (main file)
mf<- read_excel("string_interactions_short - Copy.xlsx")

## read supporting file, where matching info,is available (supporting file)
sf<- read_excel("string_mapping.xlsx")

## Get unique ids from the sf for comparison
unique_id<- sf %>%
  pull(stringId) %>% unique()

## find id that are not matched in the mapping file
unmatched_id_node1 <- mf %>% 
  filter(!(node1_string_id %in% unique_id)) %>% 
  #to remove duplicate id in case
  distinct(node1_string_id) %>%
  pull(node1_string_id)

unmatched_id_node2<- mf %>%
  filter(!(node2_string_id %in% unique_id)) %>%
  distinct(node2_string_id) %>%
  pull(node2_string_id)

## Prepare the data for excel (because of inequal data length, so adding NA for adjusting length)
max_length<- max(length(unmatched_id_node1), length(unmatched_id_node2))

## Combine these list of data into a single dataframe with two columns
unmatched_id<- data.frame(
  node1_unmatched=c(unmatched_id_node1, rep(NA, max_length-length(unmatched_id_node1))),
  node2_unmatched=c(unmatched_id_node2, rep(NA, max_length- length(unmatched_id_node2)))
)

## Save the dataframe in excel
write_xlsx(unmatched_id, path = "list_unmatched_string_id.xlsx")
