#!/bin/bash

# File containing the list of names
name_file="$1"

# File to search through
data_file="TAIR10_gene_description.txt"

# Check if the names file exists and is readable
if [[ ! -r "$name_file" ]]; then
    echo "Error: Cannot read $name_file"
    exit 1
fi

# Check if the data file exists and is readable
if [[ ! -r "$data_file" ]]; then
    echo "Error: Cannot read $data_file"
    exit 1
fi

# Loop through each name in the names file
while IFS= read -r name; do
    # Skip empty lines
    if [[ -z "$name" ]]; then
        continue
    fi
    echo "Searching for: $name" >&2
    # Use grep to search for the name in the data file
    grep -i "$name" "$data_file"
    if [ $? -eq 0 ]; then
        echo "Found: $name" >&2
    else
        echo "Not found: $name" >&2
    fi
done < "$name_file"