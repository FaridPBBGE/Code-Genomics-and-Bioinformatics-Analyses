#!/bin/bash
# Define an array of text patterns to replace
patterns=("1:chromosome_1" "2:chromosome_2" "3:chromosome_3" "4:chromosome_4" "5:chromosome_5" "6:chromosome_6" "7:chromosome_7" "8:chromosome_8" "9:chromosome_9")
# Input and output file names
input_file="Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.chr.gff3"
output_file="Zea_mays_chr_renamed.gff3"
# Loop through the patterns and run awk for each one
for pattern in "${patterns[@]}"; do
IFS=: read -ra pattern_parts <<< "$pattern"
old_text="${pattern_parts[0]}"
new_text="${pattern_parts[1]}"

awk -F'\t' -v OFS='\t' -v old="$old_text" -v new="$new_text" '{gsub(old, new, $1)}1' "$input_file" > "$output_file"

# Update the input file for the next iteration
mv "$output_file" "$input_file"
done