#!/bin/bash

mapfile -t renames < rename_sample_list.txt

# Get the count of new names
count=${#renames[@]}

# Check if the count of new names matches the count of files
if [ $count -ne $(ls -1 | grep -E '[A-Za-z0-9_-]+\.fastq\.gz' | wc -l) ]; then
  echo "Error: Number of new names does not match the number of files"
  exit 1
fi

# Iterate over the files and rename them
index=0
for file in $(ls -1 | grep -E '[A-Za-z0-9_-]+\.fastq\.gz' | sort -n); do
  oldname="$file"
  newname="${renames[$index]}"

  # Skip if the new name is empty
  if [ -z "$newname" ]; then
    echo "Skipping file $oldname as new name is empty"
    continue
  fi

  # Rename the file
  mv "$oldname" "$newname"

  echo "Renamed file $oldname to $newname"

  index=$((index + 1))
done

