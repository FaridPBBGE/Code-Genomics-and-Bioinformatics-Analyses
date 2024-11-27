#!/bin/bash
# dir containing input files
input_dir="/home/farid/Arabidopsis/gene_ID_files"
# dir for saving the output files
output_dir="/home/farid/Arabidopsis/gene_uniq_ID_files"
# loop through each file in the input dir
for file in "$input_dir"/*.txt; do
# extract base name of the file (without dir and extension)
base_name=$(basename "$file".txt)
# Defining output file name
output_file="$output_dir/${base_name}_uniq_gene.txt"
#applying the command
sort "$file" | uniq > "$output_file"
done