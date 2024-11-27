#!/bin/bash

# giving input file 
input_file="$1"
output_file="$2"

# Creating or clearing the output_file
> "$output_file"

# reading the input file one by one
while IFS= read -r var; do
echo $var
var2=$(($var + 5000))
var1=$(($var - 5000))
echo "$var:$var1;$var2" >> "$output_file"

# Process the extraction genes based on SNP position
awk -F "\t" -v start="$var1" -v end="$var2" '{ if(($1=="Chr1") && (($3<=start) && ($4>=start && $4<=end))) { print } }' TAIR10_mrna_list.txt >> "$output_file"
done < "$input_file"
