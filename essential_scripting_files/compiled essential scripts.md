## Codes for extraction information from fasta/bam/gtf/gff3 file
- **To extract genes information for certain regions of genome based on "start" and "end" column/subsetting genes list from a file for certain genomic region.**
~~~
## Below codes are for running as a script
#!/bin/bash
var=$@
echo $var
var2=$(($var+5000))
var1=$(($var-5000))
echovar "$var1;$var2";
##### based on two conditions in two columns
cat /home/farid/Vflay_genes_list | awk -F "\t" -v start="$var1" -v end="$var2" '{ if(($3>=start) && ($4<=end)) { print } }' 
##### three conditions in two columns
awk -F " " -v start="$var1" -v end="$var2" '{ if(($1=="SOVchr2") && (($3>=start && $3<=end) && ($4>=end))) { print } }' 
condition c
awk -F " " -v start="$var1" -v end="$var2" '{ if(($1=="SOVchr2") && (($3<=start) && ($4>=start && $4<=end))) { print } }' 
~~~
- **Renaming files together from a given txt file**
~~~
# Read new names from the file
mapfile -t newnames < newnames.txt
# Get the count of new names
count=${#newnames[@]}
# Check if the count of new names matches the count of files
if [ $count -ne $(ls -1 | grep -E '^[0-9]+_genes_[0-9]+kb_[abc]$' | wc -l) ]; then
  echo "Error: Number of new names does not match the number of files"
  exit 1
fi
# Iterate over the files and rename them
index=0
for file in $(ls -1 | grep -E '^[0-9]+_genes_[0-9]+kb_[abc]$' | sort -n); do
  oldname="$file"
  newname="${newnames[$index]}"
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
~~~
- **To get the read's mate unmapped from bam file** 
~~~
samtools view <file>.bam | cut –f7 | grep –c ‘*’
~~~
